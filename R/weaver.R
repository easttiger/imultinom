### Inputs
Weaver.input.validate <-
# Validates input dimensions
function(a,b,tDe){
  stopifnot("matrix" == class(tDe));
  stopifnot("numeric" == mode(tDe));
  stopifnot("numeric" == mode(a));
  stopifnot("numeric" == mode(b));
  stopifnot((ncol(tDe) == length(a)) && (nrow(tDe) == length(b)));
}

## Converter for adjacency matrix input format
Weaver.matrix <-
function(A){
  #validate A
  Weaver.input.validate.matrix(A);
  #convert A
  input = Weaver.input.matrix(A);
  #basic weaver
  WeaverBas(input$a,input$b,input$tDe)
  #superposed weaver
}

Weaver.input.validate.matrix <-
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

Weaver.input.matrix <-
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

## Special input converters
Weaver.Examples.Hunter2004PL.ConvertInput <-
function(txt="http://sites.stat.psu.edu/~dhunter/code/btmatlab/nascar2002.txt",hasHeader=TRUE){  
	dat = read.table(file=txt, header=hasHeader)
	nb = nrow(dat);
	nrace = 36; # being lazy
	ndriver = 83; # being lazy
	nsizes = tabulate(dat$Race);
	rslt = list();
	rslt$b = rep(-1,nb);  # b is always -1
	a = integer(ndriver);
	tDe = matrix(integer(nb * ndriver),ncol=ndriver);
	ib = 0;
	for(r in 1:nrace){
		nsize = nsizes[r];
		for(i in 1:nsize){
			a[dat$DriverID[ib+i]] = a[dat$DriverID[ib+i]] + 1;
			for(j in i:nsize){
				tDe[ib+i,dat$DriverID[ib+j]] = 1;
			}
		}
		ib = ib + nsize;
	}
	rslt$a = drop(a);
	rslt$tDe = tDe;
	rslt
}

Weaver.Examples.Hankin2010Volleyball.ConvertInput <-
function(n = 9){
  require("hyperdirichlet");
  data(volleyball);
  tDe = binmat(n);
  ida = tDe %*% rep(1,n) == 1L;
  b = powers(vb_synthetic);  
  a = rev(b[ida]);
  idb = !ida;
  
  idb[1] = FALSE; # the row has all 0
  idb[length(idb)] = FALSE; # the row has all 1
  idb[b == 0] = FALSE; # remove those b=0 cases
  b = b[idb];
  tDe = tDe[idb,];
  
  list(a=a,b=b,tDe=tDe)
}

Weaver <-
function(a,b,tDe, tol=1e-10,maxit=500,iteration=FALSE,ini=-1, PriorThickness=0){
  env=environment()
  res = NULL
  tryCatch({
    env$res=WeaverBas(a,b,tDe,tol=tol,maxit=maxit)
  },error = function(err){
    tryCatch({
      env$res=WeaverGre(a,b,tDe,tol=tol,maxit=maxit)
    }, error = function(err){
      env$res=WeaverBayes(a,b,tDe,tol=tol,maxit=maxit,PriorThickness=PriorThickness)
    })
  })
  if(is.null(env$res) || is.nan(env$res$iter$lnLik)){
    env$res = WeaverBayes(a,b,tDe,tol=tol,maxit=maxit,PriorThickness=PriorThickness)
  }
  env$res
}

# Algorithms
WeaverBas <-
# This function implements the Basic Weaver Algorithm
function(a,b,tDe,listinput, tol=1e-10,maxit=500,iteration=FALSE,ini=-1){
  if(!missing(listinput)){
    a = listinput$a;
    b = listinput$b;
    tDe = t(listinput$De);
  }
  m = sum(a) + sum(b);
  if(any(ini <= 0)){
    if(any(a <= 0)){
      x = rep(1,length(a))
    }else{
      x = a
    }
  }else{
    x = ini;
  }
  x = x / sum(x);
  
  iterCount = 0;
  if(iteration){
    lena = length(a);
    iter = list();
    iter$path = double(lena * maxit);
    iter$lnLik = double(maxit);
  }
  e = 1e10;
  while(e > tol && iterCount <= maxit){
    if(iteration){
      iter$path[iterCount * lena + 1:lena] = drop(x);
      iter$lnLik[iterCount+1] = sum(a * log(x)) + sum(b * log(tDe %*% x));
      
    }
    iterCount = iterCount + 1;
    
    tau = b / (tDe %*% x);
    tau0 = m - sum(tau)
    xnew = drop(a / t(tau0 + t(tau) %*% (1L - tDe)));
    
    #if(any(xnew <= 0)){
      #(-x[xnew <=0])  /  (xnew[xnew <= 0] - x[xnew <= 0]) 
    #  k = min(1 /  (1 - xnew[xnew <= 0] / x[xnew <= 0]))
    #  k = k / 3
    #  xnew = x + k * (xnew - x)
    #}
    stopifnot(all(xnew >= 0));
    #xnew[xnew <= 0] = 0.01
    
    xnew = xnew / sum(xnew);
    #e = sqrt(sum((a - x / sum(x) * drop(t(tau0 + t(tau) %*% (1L - tDe))))^2));
    e = sum(abs(x - xnew))
    #x = sapply(xnew,function(x)max(x,1e-16))
    x = xnew
    
    if(iterCount > maxit){
      warning("maxit reached");
    }
    
  }
  
  if(iteration){
    # trim 
    iter$path = matrix(iter$path[1:(lena * iterCount)], ncol=lena, byrow = TRUE);
    iter$lnLik = iter$lnLik[1:iterCount];
    iter$count = iterCount;
    
  }else{
    iter=list(
      count=iterCount,
      lnLik = sum(a * log(x)) + sum(b * log(tDe %*% x)));
  }
  list(x=x,iter=iter,e=e)
  
}

initWeaver <-
function(a,b,tDe,listinput,PriorThickness=0){
  if(!missing(listinput)){
    a = listinput$a;
    b = listinput$b;
    tDe = t(listinput$De);
  }
  if(PriorThickness == 0){
    PriorThickness = (sum(abs(b))) * (sum(abs(b))) / ((sum(abs(a))) + 1)
  }
  m1 = sum(a) + sum(b)
  m2 = m1 + PriorThickness * length(a)
  res = WeaverBayes(a,b,tDe,PriorThickness = PriorThickness,tol = 1/(length(a)^2), maxit=50,iteration=FALSE, ini=rep(1/length(a),length(a)))
  list(ini=res$x, guesssoln=(x * m2 - PriorThickness) / m1,PriThi=PriorThickness, res$iter$count)
}

WeaverBayes <-
# This function implemnets the Superposed Weaver Algorithm
function(a,b,tDe,listinput,PriorThickness=0,tol=1e-10, maxit=10000,iteration=FALSE,ini=-1){
  if(!missing(listinput)){
    a = listinput$a;
    b = listinput$b;
    tDe = t(listinput$De);
  }
  if(any(a <= 0)) maxit = max(50000,maxit)
  if(all(a == 0)){
    a = rep(1e-5,length(a))
  }
  if(any(ini <= 0)){
    if(any(a <= 0)){
      x = rep(1,length(a))
    }else{
      x = a
    }
  }else{
    x = ini;
  }
  
  x = x / sum(x);
  
  if(any(x <= 0)) x = 1.0 / length(a);
  e = 1e10;
  if(PriorThickness == 0){
    PriorThickness = (sum(abs(b))) * (sum(abs(b))) / ((sum(abs(a))) + 1) * 10
  }
  iterCount = 0;
  iterCountWeaverBayes = 0;
  if(iteration){
    iter=list();
  }
  
  while(e > tol && iterCount <= maxit){
    xnew = WeaverBas(a + x * PriorThickness,b,tDe,tol=tol,iteration=iteration, ini = a + x * PriorThickness);
    iterCount = iterCount + xnew$iter$count;
    iterCountWeaverBayes = iterCountWeaverBayes + 1;
    if(iteration){
      lnLikOffset = -sum(x * PriorThickness * log(x))
      if(length(iter) == 0){
        iter=list(path=xnew$iter$path,lnLik=xnew$iter$lnLik + lnLikOffset);
      }else{
        iter$path = rbind(iter$path, xnew$iter$path);
        iter$lnLik = c(iter$lnLik, xnew$iter$lnLik + lnLikOffset);
      }
            
    }
    xnew = xnew$x;
    stopifnot(all(xnew >= 0));
    #tau = b / (tDe %*% xnew)
    #tau0 = sum(a) + sum(b) - sum(tau)
    #e = sqrt(sum((xnew * drop(tau0 + t(tau) %*% (1L - tDe)))^2))
    e = sum(abs(x - xnew))
    
    x = sapply(xnew, function(x)max(x,1e-16))
    
    if(iterCount > maxit){
      warning("maxit reached");
    }
  }
  
  stopifnot(all(xnew >= 0));
  e = sum(abs(xnew - x));
  x = xnew;
  if(iteration){
    iter$count = list(Weaver=iterCount,WeaverBayes=iterCountWeaverBayes);
  }else{
    iter = list(
      count=list(Weaver=iterCount,WeaverBayes=iterCountWeaverBayes),
      lnLik = sum(a * log(x)) + sum(b * log(tDe %*% x)));
  }
  list(x=x,iter=iter, prithi=PriorThickness)
}

WeaverGre <-
# This function implements the Greedy Weaver Algorithm
function(a,b,tDe,listinput,tol=1e-10,maxit=500,iteration=FALSE,ini=-1){
  if(!missing(listinput)){
    a = listinput$a;
    b = listinput$b;
    tDe = t(listinput$De);
  }
  m = sum(a) + sum(b);
  if(any(ini <= 0)){
    x = WeaverBayes(a,b,tDe,tol=0.1)$x;
  }else{
    x = ini / sum(ini);
  }
  if(any(x <= 0)) x = 1.0 / length(a);
  iterCount = 0;
  if(iteration){
    lena = length(a);
    iter = list();
    iter$path = double(lena * maxit);
    iter$lnLik = double(maxit);
  }
  e = 1e10;
  dbest = 1e10;
  xbest = x;
  while(e > tol && iterCount <= maxit){
    if(iteration){
      iter$path[iterCount * lena + 1:lena] = drop(x);
      iter$lnLik[iterCount+1] = sum(a * log(x)) + sum(b * log(tDe %*% x));      
    }
    iterCount = iterCount + 1;
    
    tau = b / (tDe %*% x);
    tau0 = m - sum(tau);
    d = drop(x * (tau0 + t(1 - tDe) %*% tau) - a);
    if(sum(abs(d)) < sum(abs(dbest))){
      dbest = d;
      xbest = x;
    }
    e = sum(abs(d));
    imax = which(abs(d[1:(length(d)-1)]) == max(abs(d[1:(length(d)-1)])));
    u = c(0, x[imax], 1.05 * x[imax]);
    xtemp = x;
    xtemp[imax] = u[3];
    xtemp[length(xtemp)] = 0;
    xtemp[length(xtemp)] = 1 - sum(xtemp);
    tautemp = b / (tDe %*% xtemp);
    dtemp = drop(xtemp * (m - sum(tautemp) + t(1 - tDe) %*% tautemp) - a);
    v = c(-a[imax],d[imax],dtemp[imax]);
    gam = v[1];
    alp = (u[3] * (v[2] - v[1]) - u[2] * (v[3] - v[1])) / (u[2] * u[2] * u[3] - u[2] * u[3] * u[3]);
    bet = (-u[3] * u[3] * (v[2] - v[1]) + u[2] * u[2] * (v[3] - v[1])) / (u[2] * u[2] * u[3] - u[2] * u[3] * u[3]);
    xnew = x;
    xnew[imax] = (-bet + sqrt(bet * bet - 4.0 * alp * gam)) / (2.0 * alp);
    if(sum(xnew[1:(length(xnew)-1)]) > 1){
      xnew = xnew / sum(xnew);
    }else{
      xnew[length(xnew)] = 0;
      xnew[length(xnew)] = 1 - sum(xnew);
    }
    stopifnot(all(xnew >= 0));
    
    x = xnew;
    if(iterCount > maxit){
      warning("maxit reached");
    }
  }
  if(iteration){
    # trim 
    iter$path = matrix(iter$path[1:(lena * iterCount)], ncol=lena, byrow = TRUE);
    iter$lnLik = iter$lnLik[1:iterCount];
    iter$count = iterCount;
    
  }else{
    iter=list(
      count=iterCount,
      lnLik = sum(a * log(x)) + sum(b * log(tDe %*% x)));
  }
  
  list(x=x,iter=iter)
}

imm <-
function(p,a,b,De,listInput, logscale=TRUE){
  if(!missing(listInput)){
    a = listInput$a;
    b = listInput$b;
    De = listInput$De;
  }
  p = p / sum(p)
  loglik = sum(a * log(p)) + sum(b * log(drop(t(p) %*% De)))
  ifelse(logscale, loglik, exp(loglik))
}