setClass(Class     = "imultinom",
		 slots     = c(a="numeric",b="numeric",De="matrix"),
		 prototype = prototype(a=NA_real_, b=NA_real_, De=matrix())
		)

imultinom <-
function(listInput,matrixInput,a,b,De){
	#determine input type
	if(!missing(listInput)){
		a  = listInput$a
		b  = listInput$b
		De = listInput$De
	}else if(!missing(matrixInput)){
		.matrix.input.validate(matrixInput)
		listInput = .input.matrix(matrixInput)
		a  = listInput$a
		b  = listInput$b
		De = t(listInput$tDe)
	}

	stopifnot("matrix" == class(De))
	stopifnot("numeric" == mode(a))
	stopifnot("numeric" == mode(b))
	stopifnot(length(a) >  1L)
	stopifnot(length(b) >= 1L)
	stopifnot((nrow(De) == length(a)) && (ncol(De) == length(b)))
	stopifnot(all(a >= 0))

	new("imultinom",a=a,b=b,De=De)
}

setGeneric("likelihood",function(x, p,logscale=FALSE){standardGeneric("likelihood")})
setMethod("likelihood", signature="imultinom",
function(x, p,logscale){
  p = p / sum(p)
  loglik = sum(x@a * log(p)) + sum(x@b * log(drop(t(p) %*% x@De)))
  ifelse(logscale, loglik, exp(loglik))
})

setGeneric("mle",function(x){standardGeneric("mle")})
setMethod("mle", signature="imultinom",
function(x){
	Weaver(x@a,x@b,t(x@De),tol = getOption("tolerance"), iteration = getOption("iterrec"),maxit = getOption("maxit"),PriorThickness = getOption("PriorThickness"))
})

".matrix.input.validate" <-
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

".input.matrix" <-
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

