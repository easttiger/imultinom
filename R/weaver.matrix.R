## Converter for adjacency matrix input format
weaver.matrix <-
function(A){
  #validate A
  weaver.matrix.input.validate(A);
  #convert A
  input = weaver.input.matrix(A);
  #basic weaver
  WeaverBas(input$a,input$b,input$tDe)
  #superposed weaver
}

weaver.matrix.input.validate <-
# validate certain features in the adjacency matrix input
function(A){
  nr = nrow(A);
  nc = ncol(A);
  #stops
  stopifnot(nr >= 2);
  stopifnot(nr == nc);
  stopifnot(all(diag(A) == 0));
  
  #warnings
  if(any(A < 0)){
    warning('Negative element detected in the adj. mat. input')
  }  
  #potentially validate strong connectivity
}

weaver.input.matrix <-
# converts adjacency matrix input to list(a,b,tDe)
function(A){
  n = nrow(A);
  rslt = list();
  rslt$a = drop(A %*% rep(1, n));
  b = double(n * (n - 1) / 2);
  De = integer(n * n * (n - 1) / 2);
  nb = 0;
  for(i in seq(2,n)){
    for(j in seq(1,i-1)){
      s = A[j,i] + A[i,j];
      if( s > 0 ){        
        b[nb + 1] = -s;
        De[nb * n + i] = 1L;
        De[nb * n + j] = 1L;
        nb = nb + 1;
      }
    }
  }
  if(nb == 0){
    rslt$b = 0;
    rslt$tDe = t(rep(0L,n));
  }else{
    rslt$b = b[1:nb];
    rslt$tDe = matrix(De[1:(n * nb)], 
              ncol = n, byrow = TRUE);
  }
  rslt
}
