Weaver <-
function(a,b,tDe, tol=1e-10,maxit=500,iteration=FALSE,ini=-1, PriorThickness=0){
  Weaver.input.validate(a,b,tDe)
  env=environment()
  res = NULL
  if(!is.logical(iteration)) iteration = FALSE
  if(is.null(maxit)){maxit = 500}
  if(is.null(PriorThickness)){PriorThickness = 0}
  if(is.null(tol)){tol = 1e-10}
  tryCatch({
    env$res=weaver.vanilla(a,b,tDe,tol=tol,maxit=maxit,iteration = iteration)
  },error = function(err){
    tryCatch({
      env$res=weaver.greedy(a,b,tDe,tol=tol,maxit=maxit,iteration = iteration)
    }, error = function(err){
      env$res=weaver.bayes(a,b,tDe,tol=tol,maxit=maxit,PriorThickness=PriorThickness,iteration = iteration)
    })
  })
  if(is.null(env$res) || is.nan(env$res$iter$lnLik)){
    env$res = weaver.bayes(a,b,tDe,tol=tol,maxit=maxit,PriorThickness=PriorThickness,iteration = iteration)
  }
  env$res
}

# Algorithms
weaver.vanilla <-
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

    stopifnot(all(xnew >= 0));


    xnew = xnew / sum(xnew);
    if(iteration){
      sse = sqrt(sum((a - x / sum(x) * drop(t(tau0 + t(tau) %*% (1L - tDe))))^2));
      iter$sse = c(iter$sse, sse)
    }

    e = sum(abs(x - xnew)) #another error function
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

  list(x=x,iter=iter,e=e,covmle=cov.of.mle(x,a,b,t(tDe)))

}

weaver.bayes <-
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
  iterCountweaver.bayes = 0;
  if(iteration){
    iter=list();
  }

  while(e > tol && iterCount <= maxit){
    xnew = weaver.vanilla(a + x * PriorThickness,b,tDe,tol=tol,iteration=iteration, ini = a + x * PriorThickness);
    iterCount = iterCount + xnew$iter$count;
    iterCountweaver.bayes = iterCountweaver.bayes + 1;
    if(iteration){
      lnLikOffset = -sum(x * PriorThickness * log(x))
      if(length(iter) == 0){
        iter=list(path=xnew$x, fullpath=xnew$iter$path,lnLik=xnew$iter$lnLik[length(xnew$iter$lnLik)]+ lnLikOffset, fulllnLik=xnew$iter$lnLik + lnLikOffset, sse=sse(xnew$x, a, b, tDe));
      }else{
        iter$path = rbind(iter$path,xnew$x)
        iter$fullpath = rbind(iter$fullpath, xnew$iter$path);
        iter$fulllnLik = c(iter$fulllnLik, xnew$iter$lnLik + lnLikOffset);
        iter$lnLik = c(iter$lnLik, iter$fulllnLik[length(iter$fulllnLik)])
        iter$sse = c(iter$sse, sse(xnew$x, a, b, tDe))
      }

    }
    xnew = xnew$x;
    stopifnot(all(xnew >= 0));
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
    iter$count = list(Weaver=iterCount,weaver.bayes=iterCountweaver.bayes);
  }else{
    iter = list(
      count=list(Weaver=iterCount,weaver.bayes=iterCountweaver.bayes),
      lnLik = sum(a * log(x)) + sum(b * log(tDe %*% x)));
  }
  list(x=x,iter=iter, e=e, prithi=PriorThickness,covmle=cov.of.mle(x,a,b,t(tDe)))
}

weaver.greedy <-
# This function implements the Greedy Weaver Algorithm
function(a,b,tDe,listinput,tol=1e-10,maxit=500,iteration=FALSE,ini=-1){
  if(!missing(listinput)){
    a = listinput$a;
    b = listinput$b;
    tDe = t(listinput$De);
  }
  m = sum(a) + sum(b);
  if(any(ini <= 0)){
    x = weaver.bayes(a,b,tDe,tol=0.1)$x;
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
    if(iteration){
      sse = sqrt(sum((a - x / sum(x) * drop(t(tau0 + t(tau) %*% (1L - tDe))))^2));
      iter$sse = c(iter$sse, sse)
    }
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

  list(x=x,iter=iter,e=e,covmle=cov.of.mle(x,a,b,t(tDe)))
}

## Interface functions
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
  res = weaver.bayes(a,b,tDe,PriorThickness = PriorThickness,tol = 1/(length(a)^2), maxit=50,iteration=FALSE, ini=rep(1/length(a),length(a)))
  list(ini=res$x, guesssoln=(x * m2 - PriorThickness) / m1,PriThi=PriorThickness, res$iter$count)
}

Weaver.input.validate <-
# Validates input dimensions
function(a,b,tDe){
  stopifnot("matrix" == class(tDe));
  stopifnot("numeric" == mode(a));
  stopifnot("numeric" == mode(b));
  stopifnot((ncol(tDe) == length(a)) && (nrow(tDe) == length(b)));
}


obs.info.mat <-
function(p,a,b,De){
	d = length(p)
	q = length(b)
	J = cbind(diag(rep(1,d-1)),rep(-1,d-1))
	dDe = J %*% De
	tDep2 = as.vector((t(De) %*% as.vector(p)))^2
	psi = matrix(double((d-1)^2),nrow=d-1)
	for(i in 1:(d-1)){
	  for(k in i:(d-1)){
		psi[i,k] = sum(b * dDe[i,] * dDe[k,] / tDep2)
		if(i != k) psi[k,i] = psi[i,k]
	  }
	}
	diag((a / (p^2))[-d]) + a[d] / (p[d]^2) + psi
}

cov.of.mle <-
function(p,a,b,De){
	covExLast <- solve(obs.info.mat(p,a,b,De))
	covOfLast <- as.vector(-apply(covExLast,1,sum))
	varOfLast <- sum(covExLast)
	cov.of.mle <- rbind(cbind(covExLast,covOfLast),cbind(t(covOfLast),varOfLast))
	dimnames(cov.of.mle) <- NULL
	cov.of.mle
}

sse <-
function(x,a,b,tDe){
  tau = b / (tDe %*% x);
  tau0 = sum(a) + sum(b) - sum(tau)
  sqrt(sum((a - x / sum(x) * drop(t(tau0 + t(tau) %*% (1L - tDe))))^2))
}
