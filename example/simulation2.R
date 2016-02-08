# setwd("E:/Dropbox/Subjects/W/Weaver Paper Reading list/Real Data Examples")
source("weaver.R")
RUN_PARALLEL = T
library(snow)
m=1;
randparti <-
  function(p){
    stopifnot(all(p>0))
    p = p / sum(p)
    d = length(p)
    icol = 0;
    #draw nu from U{10d, 10d+1, ..., 100d}
    #m_range = seq(from= 10 * d, to= 100 * d, by= 1L)
    m_range = seq(from= m * d, to= 2 * m * d, by= 1L)
    nu = sample(x = m_range, size = 1L)
    
    #draw partition
    #draw an upper bound of the number of subsets in the partition
    nSubsets = sample(x = 2:d, size=1L)
    #draw subset membership for atoms 1..d
    #it may not be that all subsets turned up
    #include 0 for the simulation of conditioning
    mSubsets = sample(x = 0:nSubsets, size=d, replace=T)
    uSubsets = unique(sort(mSubsets))
    nSubsets = length(uSubsets) # if 0 is included in the outcome, 
    # there will be an additional negative-exponent subset, 
    # but that will not break this expression
    
    #if(RUN_PARALLEL){
    #cl <- makeSOCKcluster(rep("localhost",ncluster))
    #clusterExport(cl, list("p","mSubsets","d"))
    #p_Xi = parSapply(cl, uSubsets,function(x){ifelse(x==0,0,sum(p[which(mSubsets == x)]))})
    
    #}else{
    p_Xi = sapply(uSubsets,function(x){ifelse(x==0,0,sum(p[which(mSubsets == x)]))})
    #}
    
    #if(RUN_PARALLEL){
    #  Xi = parSapply(cl,uSubsets,function(x){
    #    res=rep(0,d);
    #    if(x>0){
    #      res[which(mSubsets == x)] = 1
    #    };
    #    res})
    #stopCluster(cl)
    #}else{
    Xi = sapply(uSubsets,function(x){
      res=rep(0,d);
      if(x>0){
        res[which(mSubsets == x)] = 1
      };
      res})
    #}
    y = as.vector(rmultinom(1,nu,p_Xi))
    
    if(uSubsets[1] == 0){
      #need to modify the first col of Xi and first element of b
      Xi[,1] = apply(Xi,1,sum)
      y[1] = -sum(y)
    }
    list(y=y,Xi=Xi,nu=nu)
  }

g <-  
  function(p,R,includeCompleteData = TRUE){
    #if(RUN_PARALLEL) cl <- makeSOCKcluster(c("localhost","localhost"))
    p = p / sum(p)
    d = length(p)
    stopifnot(R>=1L)
    #if(RUN_PARALLEL){
    #cl <- makeSOCKcluster(rep("localhost",2))
    #clusterExport(cl, list("p","d","randparti","RUN_PARALLEL","cl"))
    #res = parLapply(cl, 1:R, function(x)randparti(p))
    #stopCluster(cl)
    #}else{
    res = lapply(1:R, function(x)randparti(p))
    #}
    y = Reduce("c", sapply(res, function(x)x$y))
    Xi = Reduce("cbind", sapply(res, function(x)x$Xi))
    nu = Reduce("+", sapply(res, function(x)x$nu))
    Xi = matrix(Xi, nrow = d)
    #index of a
    colsum = as.vector(apply(Xi,2,sum))
    
    idx_b = which((colsum > 1) & (colsum < d))
    b = y[idx_b]
    De = Xi[,idx_b]
    
    if(includeCompleteData){
      a = as.vector(apply(rmultinom(R, nu, p), 1, sum)) #initialise a with prior thickness
    }else{
      a = double(d) + 0.00001
    }
    if(length(y[colsum==1]) > 1){
      a = a + as.vector(Xi[,colsum == 1] %*% y[colsum == 1]); 
    }
    if(length(y[colsum==1]) == 1){
      a = a + as.vector(Xi[,colsum == 1] * y[colsum == 1]); 
    }
    #if(RUN_PARALLEL) stopCluster(cl)
    list(a=a,b=b,De=De)
    #print(c(nu,sum(a),sum(b[b<0]),sum(b[b>0])))
  }

mleSim <-
  function(samplesize, p, R, ncluster=1, incComplete=T){
    #if(RUN_PARALLEL) cl <- makeSOCKcluster(c("localhost","localhost"))
    stopifnot(all(p>0))
    d = length(p)
    p = p / sum(p)
    
    #  res = matrix(double(d * samplesize), nrow=d)
    #  for(i in 1:samplesize){
    #    lik=g(p,R);
    #    res[,i] = as.vector(WeaverBayes(lik$a,lik$b,t(lik$De),PriorThickness = prithick, maxit = maxiter)$x);
    #  }
    if(ncluster > 1){
      require(snow)
      #stopCluster(cl)
      cl <- makeSOCKcluster(rep("localhost",ncluster))
      clusterExport(cl, list("p","d","R"), envir = sys.frame(sys.nframe()))
      clusterExport(cl, list("randparti","incComplete","Weaver","WeaverBayes","WeaverBas","WeaverGre","initWeaver","g","m"), envir = .GlobalEnv)
      
      phat = parSapply(cl, 1:samplesize, function(x){
        lik=g(p,R,includeCompleteData = incComplete);
        if(m * 1.7^R >=d){
          WeaverBas(lik$a,lik$b,t(lik$De))$x
        }else{
          Weaver(lik$a,lik$b,t(lik$De))$x
        }
      })
      
      stopCluster(cl)
    }else{
      phat = sapply(1:samplesize, function(x){
        lik=g(p,R, includeCompleteData = incComplete);
        if(m>=d/2){
          WeaverBas(lik$a,lik$b,t(lik$De))$x
        }else{
          Weaver(lik$a,lik$b,t(lik$De))$x
        }
      })
    }
    #if(RUN_PARALLEL) stopCluster(cl)
    #oldmfrow=par("mfrow")
    #par(mfrow=c(3,1))
    #plot(p,apply(phat,1,mean));abline(0,1)
    #boxplot(t(phat));
    #plot(p,apply(phat,1,var));
    #par(mfrow=oldmfrow)
    
    phat
  }

function(){
  p = 1:100; p = p/sum(p);d=length(p)
  lastvar=matrix(double(600*10),ncol=600);
  for(i in 1:10){
    phat<-mleSim(600,p,i*2,12);
    lastvar[i,] = phat[d,]
  }
  par(mfrow=c(1,1))
  par(mar=c(4,5.1,1.5,1.5))
  par(mgp=c(3.5,0.6,0))
  #boxplot(t(lastvar),names=seq(2,20,2),main=expression(bold(paste("Decreasing var(",hat(italic(p))[100],") as ", italic(R)," increases",sep=""))),ylab=expression(hat(italic(p))[100]))
  boxplot(t(lastvar),names=seq(2,20,2),main="",ylab=expression(hat(italic(p))[100]))
  title(xlab=expression(bold(R)),line=2)
  
  phatmat = matrix(double(4*100),nrow=100)
  i=1
  for(R in c(2,4,8,20)){
    1->m;
    phatmat[,i] <- mleSim(2,p,R)[,1]
    i = i + 1;
  }
  
  par(mfrow=c(2,2))
  i = 1
  for(R in c(2,4,8,20)){
    plot(p, phatmat[,i], main=paste("m = ",m, ", R = ",R, sep=""), ylab=expression(hat(bolditalic(p))),xlab=""); abline(0,1)
    title(xlab=expression(bolditalic(p)),line=2)
    i = i + 1
  }
  
}

function(){
  
  p=1:100;p=p/sum(p);
  par(mfrow=c(2,2));
  par(mar=c(4,5.1,1.5,1.5))
  par(mgp=c(3.5,0.6,0))
  for(i in c(1,10,100,1000)){
    i->m;
    phat=mleSim(2,p,2,2);
    plot(p,phat[,1],main=paste("m = ",i, ", R = ",2, sep=""), ylab=expression(hat(bolditalic(p))),xlab="");abline(0,1)
    title(xlab=expression(bolditalic(p)),line=2)
  }
  
}

function(){
  # 4 x 4
  p=1:100;p=p/sum(p);
  par(mfrow=c(4,4));
  par(mar=c(4,5.1,1.5,1.5))
  par(mgp=c(3.5,0.6,0))
  for(m in c(1,10,100,1000)){
    for(R in c(2,4,8,16)){
      phat=mleSim(2,p,R,2);
      plot(p,phat[,1],main=paste("m = ",m, ", R = ",R, sep=""), ylab=expression(hat(bolditalic(p))),xlab="");abline(0,1)
      title(xlab=expression(bolditalic(p)),line=2)
    }
  }
}

function(){
  p = 1:100; p = p/sum(p);d=length(p)
  2->m;
  lastvar=matrix(double(600*10),ncol=600);
  for(i in 1:10){
    phat<-mleSim(600,p,i*2,12);
    lastvar[i,] = phat[d,]
  }
  par(mfrow=c(1,1))
  par(mar=c(4,5.1,1.5,1.5))
  par(mgp=c(3.5,0.6,0))
  #boxplot(t(lastvar),names=seq(2,20,2),main=expression(bold(paste("Decreasing var(",hat(italic(p))[100],") as ", italic(R)," increases",sep=""))),ylab=expression(hat(italic(p))[100]))
  boxplot(t(lastvar),names=seq(2,20,2),main=paste("m = ",m, ", sample size = 600", sep=""),ylab=expression(hat(italic(p))[100]))
  title(xlab=expression(bold(R)),line=2)
}

function(){
  p = 1:100; p = p/sum(p);d=length(p)
  R = 2;
  
  lastvar=matrix(double(600*10),ncol=600);
  for(i in 1:10){
    (2^i)->m;
    phat<-mleSim(600,p,R,12);
    lastvar[i,] = phat[d,]
  }
  par(mfrow=c(1,1))
  par(mar=c(4,5.1,1.5,1.5))
  par(mgp=c(3.5,0.6,0))
  #boxplot(t(lastvar),names=seq(2,20,2),main=expression(bold(paste("Decreasing var(",hat(italic(p))[100],") as ", italic(R)," increases",sep=""))),ylab=expression(hat(italic(p))[100]))
  boxplot(t(lastvar),names=(2^(1:10)),main=paste("R = ",R, ", sample size = 600", sep=""),ylab=expression(hat(italic(p))[100]))
  title(xlab=expression(bold(m)),line=2)
}

function(){
  p = 1:100; p = p/sum(p);d=length(p)
  2->m;
  firstvar=matrix(double(600*10),ncol=600);
  for(i in 1:10){
    phat<-mleSim(600,p,i*2,12);
    firstvar[i,] = phat[1,]
  }
  par(mfrow=c(1,1))
  par(mar=c(4,5.1,1.5,1.5))
  par(mgp=c(3.5,0.6,0))
  #boxplot(t(lastvar),names=seq(2,20,2),main=expression(bold(paste("Decreasing var(",hat(italic(p))[100],") as ", italic(R)," increases",sep=""))),ylab=expression(hat(italic(p))[100]))
  boxplot(t(firstvar),names=seq(2,20,2),main=paste("m = ",m, ", sample size = 600", sep=""),ylab=expression(hat(italic(p))[1]))
  title(xlab=expression(bold(R)),line=2)
}

function(){
  p = 1:100; p = p/sum(p);d=length(p)
  R = 2;
  
  firstvar=matrix(double(600*10),ncol=600);
  for(i in 1:10){
    (2^i)->m;
    phat<-mleSim(600,p,R,12);
    firstvar[i,] = phat[1,]
  }
  par(mfrow=c(1,1))
  par(mar=c(4,5.1,1.5,1.5))
  par(mgp=c(3.5,0.6,0))
  #boxplot(t(lastvar),names=seq(2,20,2),main=expression(bold(paste("Decreasing var(",hat(italic(p))[100],") as ", italic(R)," increases",sep=""))),ylab=expression(hat(italic(p))[100]))
  boxplot(t(firstvar),names=(2^(1:10)),main=paste("R = ",R, ", sample size = 600", sep=""),ylab=expression(hat(italic(p))[1]))
  title(xlab=expression(bold(m)),line=2)
}

function(){
  #2x2 varplot
  p = 1:100; p = p/sum(p);d=length(p)
  R = 2;
  par(mfcol=c(2,2))
  par(mar=c(4,5.1,1.5,1.5))
  par(mgp=c(3.5,0.6,0))
  
  firstvar=matrix(double(600*10),ncol=600);
  for(i in 1:10){
    (2^i)->m;
    phat<-mleSim(600,p,R,12);
    firstvar[i,] = phat[1,]
    lastvar[i,] = phat[d,]
  }
  
  #boxplot(t(lastvar),names=seq(2,20,2),main=expression(bold(paste("Decreasing var(",hat(italic(p))[100],") as ", italic(R)," increases",sep=""))),ylab=expression(hat(italic(p))[100]))
  boxplot(t(firstvar),names=(2^(1:10)),main=paste("R =",R),ylab=expression(hat(italic(p))[1]))
  title(xlab=expression(bold(m)),line=2)
  
  #boxplot(t(lastvar),names=seq(2,20,2),main=expression(bold(paste("Decreasing var(",hat(italic(p))[100],") as ", italic(R)," increases",sep=""))),ylab=expression(hat(italic(p))[100]))
  boxplot(t(lastvar),names=(2^(1:10)),main=paste("R =",R),ylab=expression(hat(italic(p))[100]))
  title(xlab=expression(bold(m)),line=2)
  
  
  2->m;
  firstvar=matrix(double(600*10),ncol=600);
  for(i in 1:10){
    phat<-mleSim(600,p,i*2,12);
    firstvar[i,] = phat[1,]
    lastvar[i,] = phat[d,]
  }
  #boxplot(t(lastvar),names=seq(2,20,2),main=expression(bold(paste("Decreasing var(",hat(italic(p))[100],") as ", italic(R)," increases",sep=""))),ylab=expression(hat(italic(p))[100]))
  boxplot(t(firstvar),names=seq(2,20,2),main=paste("m =",m),ylab=expression(hat(italic(p))[1]))
  title(xlab=expression(bold(R)),line=2)
  #boxplot(t(lastvar),names=seq(2,20,2),main=expression(bold(paste("Decreasing var(",hat(italic(p))[100],") as ", italic(R)," increases",sep=""))),ylab=expression(hat(italic(p))[100]))
  boxplot(t(lastvar),names=seq(2,20,2),main=paste("m =",m),ylab=expression(hat(italic(p))[100]))
  title(xlab=expression(bold(R)),line=2)
  
}

function(){
  # 4 x 4
  p=1:100;p=p/sum(p);
  par(mfrow=c(2,2));
  par(mar=c(4,5.1,1.5,1.5))
  par(mgp=c(3.5,0.6,0))
  for(m in c(100,1000)){
    for(R in c(100,1000)){
      phat=mleSim(2,p,R,2);
      plot(p,phat[,1],main=paste("m = ",m, ", R = ",R, sep=""), ylab=expression(hat(bolditalic(p))),xlab="");abline(0,1)
      title(xlab=expression(bolditalic(p)),line=2)
    }
  }
  
}
require(ggplot2)
function(){
  p = c(1, rep(100,99)); p=p/sum(p);
  a = c(1, rep(100,99))
  a[1] = a[1] + 1
  b = c(101, rep(200,49), c(rep(200,48),300))
  De = matrix(double(100*99),nrow =100); 
  for(i in 1:50){s = (i-1)*2+1; De[c(s,s+1),i]=1}
  for(i in 51:99){s = (i-51)*2+2; De[c(s,s+1),i]=1}
  De[100,99] = 1
  WeaverBayes(a,b,t(De),tol=1e-15, maxit=10000, iteration=T)->res
  
  layout(matrix(c(1,1,1,2,2,3), 2, 3, byrow = TRUE))
  par(mar=c(2.4,2.7,2.4,1.5))
  par(mgp=c(1.5,0.4,0))
  plot(res$iter$lnLik,xlab="Iterations",ylab="Log-likelihood", xaxt='n', type='p', cex=0.3)
  axis(1, at=c(1,seq(500,6500,1000),7038))
  
  plot(res$iter$lnLik[1:250],xlab="Iterations",ylab="Log-likelihood", type='p', cex=0.8, xaxt='n', ylim = c(-122724,-122723.4))
  axis(1, at=c(1,seq(50,250,50)))
  
  plot(res$iter$lnLik[7000:7038],xlab="Iterations",ylab="Log-likelihood", type='p', cex=0.9, xaxt='n',yaxt='n', ylim = c(-122724,-122723.4))
  axis(1, at=c(0,seq(10,30,10),38), labels=c(7000,seq(7010,7030,10),7038))
  axis(2, at=c(-122723.8, -122723.6, -122723.4))
  
}

function(){
  #this function makes a plot of the simple example 2D
  a = c(2,2,2)
  b = 4
  De = c(1,1,0)
  it = WeaverBas(a,b,t(De),iteration=T)
  par(mfrow=c(1,3))
  
  plot(it$iter$path[1:5,1]*10,it$iter$path[1:5,2]*10,type='p',cex=1,pch=19, xlab=expression(x[1]),ylab=expression(x[2]))
  text(it$iter$path[1:5,1]*10,it$iter$path[1:5,2]*10, label=0:4,pos=c(4,1,2,3,4))
  points(4,4,pch=10,col="red",cex=3)
  
  plot(it$iter$path[1:5,1]*10,it$iter$path[1:5,3]*10,type='p',cex=1,pch=19, xlab=expression(x[1]),ylab=expression(x[3]))
  text(it$iter$path[1:5,1]*10,it$iter$path[1:5,3]*10, label=0:4,pos=c(1,2,3,4,1))
  points(4,2,pch=10,col="red",cex=3)
  
  plot(it$iter$path[1:5,2]*10,it$iter$path[1:5,3]*10,type='p',cex=1,pch=19, xlab=expression(x[2]),ylab=expression(x[3]))
  text(it$iter$path[1:5,2]*10,it$iter$path[1:5,3]*10, label=0:4,pos=c(1,2,3,4,1))
  points(4,2,pch=10,col="red",cex=3)
  
}


#concrete example fixed partition

sim_concreteExampleFixedPartition <-
  function(samplesize=1){
    # this block simulates from 1.2 A concrete example
    # fix: partition and sample sizes
    res = NULL
    chi2 = NULL
    phat = NULL
    S = NULL
    sdphat = NULL
    covar = NULL
    covphat = NULL
    p = t(matrix(c(0.1654,0.2024,0.1444,0.1532,0.2301,0.1046),nrow=3))
    De = matrix(c(1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1),nrow=6)
    n = c(120,40,40,100)
    for(i in 1:samplesize){
      y = list()
      y[[1]] = as.vector(table(sample(x=c(1,4,2,5,3,6),prob=p,size=n[1],replace=T)))
      y[[2]] = as.vector(table(sample(x=c("135","246"),prob=apply(p,1,sum),size=n[2],replace=T)))
      y[[3]] = as.vector(table(sample(x=c("14","25","36"),prob=apply(p,2,sum),size=n[3],replace=T)))
      y[[3]][1] = y[[3]][1] - 100
      y[[4]] = as.vector(table(sample(x=c("1","4"),prob=p[,1],size=n[4],replace=T)))
      a = y[[1]]
      a[1] = a[1] + y[[4]][1]
      a[4] = a[4] + y[[4]][2]
      b = c(y[[2]], y[[3]])
      phat = Weaver(a,b,t(De))$x
      obsInfoMat = oim(phat, a,b,De)
      sd = sqrt(as.vector(diag(covMLE(phat,a,b,De))))
      covar = covMLE(phat,a,b,De)
      if(is.null(res)){
        res = phat
        chi2 = as.double(t(t(p)[-6] - phat[-6]) %*% obsInfoMat %*% (t(p)[-6] - phat[-6]))
        sdphat = sd
        covphat = covar
      }else{
        res = c(res,phat)
        chi2 = c(chi2,as.double(t(t(p)[-6] - phat[-6]) %*% obsInfoMat %*% (t(p)[-6] - phat[-6])))
        sdphat = c(sdphat,sd)
        covphat = covphat + covar
      }
    }
    list(phat=matrix(res,nrow=6),sdphat=matrix(sdphat,nrow=6),chi2=as.vector(chi2),sephat = covphat/samplesize)
  }


oim <-
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

covMLE <-
  function(p,a,b,De){
    covExLast <- solve(oim(p,a,b,De))
    covOfLast <- as.vector(-apply(covExLast,1,sum))
    varOfLast <- sum(covExLast)
    covMLE <- rbind(cbind(covExLast,covOfLast),cbind(t(covOfLast),varOfLast))
    dimnames(covMLE) <- NULL
    covMLE
  }

coverage <-
  function(n){
    y=sim_concreteExampleFixedPartition(n)
    z=apply(y$phat,2,function(x)(x-p))/y$sdphat
    proba = seq(0.01,0.99,0.01)
    plot(proba,sapply(proba, function(x){sum(y$chi2 < qchisq(x,5)) / n}),main="Coverage Probability Plot", xlab="Prespecified Coverage Probability", ylab="Sampling Coverage Probability")
    abline(0,1)
  }


sim_concreteExampleFixedPartitionNoCompleteDataUnindentifiable <-
  function(n){
    res = NULL
    chi2 = NULL
    phat = NULL
    S = NULL
    sdphat = NULL
    covar = NULL
    covphat = NULL
    
    p = c(0.1654,0.2024,0.1444,0.1532,0.2301,0.1046)
    #p = (1:6) / sum(1:6)
    #sample2 - marginal gender
    p2 = c(sum(p[1:3]),sum(p[4:6]))
    ssize2 = 40
    
    #sample3 - marginal age
    p3 = c(p[1]+p[4],p[2]+p[5],p[3]+p[6])
    ssize3 = 40
    
    #sample4 - conditional
    p4 = c(p[1],p[4])
    ssize4 = 100
    
    
    #sample5 - conditional on male
    p5 = p[4:6]
    ssize5 = 40
    
    De = matrix(c(1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1),nrow=6)
    
    cl <- makeSOCKcluster(rep("localhost",12))
    clusterExport(cl, list("p","p2","p3","p4","p5","ssize2","ssize3","ssize4","ssize5","De","Weaver","WeaverBas","WeaverBayes","WeaverGre","cl","oim","covMLE"),envir=environment())
    res = parLapply(cl, 1:n, function(x){
      x2 = as.vector(table(sample(x=c("123","456"),prob = p2, replace = T, size = ssize2)))
      x3 = as.vector(table(sample(x=c("14","25","36"),prob = p3, replace = T, size = ssize3)))
      x4 = as.vector(table(sample(x=c("1","4"),prob = p4, replace = T, size = ssize4)))
      x5 = as.vector(table(sample(x=c("4","5","6"),prob = p5, replace = T, size = ssize5)))
      x5 = c(0,0,0)
      a = c(x4[1],0,0,x4[2]+x5[1],x5[2],x5[3])
      b = c(x2[1],x2[2]-sum(x5),x3[1]-sum(x4),x3[2],x3[3])
      phat = Weaver(a,b,t(De))$x
      a = a + 1e-10
      obsInfoMat = oim(phat, a,b,De)
      sd = sqrt(as.vector(diag(covMLE(phat,a,b,De))))
      covar = covMLE(phat,a,b,De)
      chi2 = as.double(t(t(p)[-6] - phat[-6]) %*% obsInfoMat %*% (t(p)[-6] - phat[-6]))
      list(phat=phat,sd=sd,covar=covar,chi2 = chi2)
    })
    stopCluster(cl)
    
    res
  }

sim_concreteExampleFixedPartitionNoCompleteDataModif1 <-
  function(n){
    res = NULL
    chi2 = NULL
    phat = NULL
    S = NULL
    sdphat = NULL
    covar = NULL
    covphat = NULL
    
    p = c(0.1654,0.2024,0.1444,0.1532,0.2301,0.1046)
    #p = (1:6) / sum(1:6)
    #sample2 - marginal gender
    p2 = c(sum(p[1:3]),sum(p[4:6]))
    ssize2 = 40
    
    #sample3 - marginal age
    p3 = c(p[1]+p[4],p[2]+p[5],p[3]+p[6])
    ssize3 = 40
    
    #sample4 - conditional
    p4 = c(p[1],p[4])
    ssize4 = 100
    
    
    #sample5 - conditional on male
    p5 = c(p[2]+p[6], p[3]+p[5])
    ssize5 = 40
    print(environment())
    De = matrix(c(c(1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1),c(0,1,0,0,0,1,0,0,1,0,1,0,0,1,1,0,1,1)),nrow=6)
    
    cl <- makeSOCKcluster(rep("localhost",12))
    clusterExport(cl, list("p","p2","p3","p4","p5","ssize2","ssize3","ssize4","ssize5","De","Weaver","WeaverBas","WeaverBayes","WeaverGre","cl","oim","covMLE"),envir=environment())
    res = parLapply(cl, 1:n, function(x){
      x2 = as.vector(table(sample(x=c("123","456"),prob = p2, replace = T, size = ssize2)))
      x3 = as.vector(table(sample(x=c("14","25","36"),prob = p3, replace = T, size = ssize3)))
      x4 = as.vector(table(sample(x=c("1","4"),prob = p4, replace = T, size = ssize4)))
      x5 = as.vector(table(sample(x=c("26","35"),prob = p5, replace = T, size = ssize5)))
      a = c(x4[1],0,0,x4[2],0,0)
      b = c(x2[1],x2[2],x3[1]-sum(x4),x3[2],x3[3],x5,-sum(x5))
      phat = Weaver(a,b,t(De))$x
      obsInfoMat = oim(phat, a,b,De)
      sd = sqrt(as.vector(diag(covMLE(phat,a,b,De))))
      covar = covMLE(phat,a,b,De)
      chi2 = as.double(t(t(p)[-6] - phat[-6]) %*% obsInfoMat %*% (t(p)[-6] - phat[-6]))
      list(phat=phat,sd=sd,covar=covar,chi2 = chi2)
    })
    stopCluster(cl)
    
    res
  }


sim_concreteExampleFixedPartitionNoCompleteDataModif2 <-
  function(n){
    res = NULL
    chi2 = NULL
    phat = NULL
    S = NULL
    sdphat = NULL
    covar = NULL
    covphat = NULL
    
    p = c(0.1654,0.2024,0.1444,0.1532,0.2301,0.1046)
    #p = (1:6) / sum(1:6)
    #sample2 - marginal gender
    p2 = c(sum(p[1:3]),sum(p[4:6]))
    ssize2 = 40
    
    #sample3 - marginal age
    p3 = c(p[1]+p[4],p[2]+p[5],p[3]+p[6])
    ssize3 = 40
    
    #sample4 - conditional
    p4 = c(p[1],p[4])
    ssize4 = 100
    
    
    #sample5 - conditional on male
    p5 = c(p[2]+p[6], p[3]+p[5])
    ssize5 = 40
    
    #sample6 - replace sample4
    p6 = c(p[1]+p[5], p[2]+p[4])
    ssize6 = 40
    
    print(environment())
    De = matrix(c(c(1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1),c(0,1,0,0,0,1,0,0,1,0,1,0,0,1,1,0,1,1),c(1,0,0,0,1,0,0,1,0,1,0,0,1,1,0,1,1,0)),nrow=6)
    
    cl <- makeSOCKcluster(rep("localhost",12))
    clusterExport(cl, list("p","p2","p3","p4","p5","p6","ssize2","ssize3","ssize4","ssize5","ssize6","De","Weaver","WeaverBas","WeaverBayes","WeaverGre","cl","oim","covMLE"),envir=environment())
    res = parLapply(cl, 1:n, function(x){
      x2 = as.vector(table(sample(x=c("123","456"),prob = p2, replace = T, size = ssize2)))
      x3 = as.vector(table(sample(x=c("14","25","36"),prob = p3, replace = T, size = ssize3)))
      
      x5 = as.vector(table(sample(x=c("26","35"),prob = p5, replace = T, size = ssize5)))
      x6 = as.vector(table(sample(x=c("15","24"),prob = p6, replace = T, size = ssize6)))
      a = c(0,0,0,0,0,0)
      b = c(x2[1],x2[2],x3[1],x3[2],x3[3],x5,-sum(x5),x6,-sum(x6))
      phat = Weaver(a,b,t(De))$x
      obsInfoMat = oim(phat, a,b,De)
      sd = sqrt(as.vector(diag(covMLE(phat,a,b,De))))
      covar = covMLE(phat,a,b,De)
      chi2 = as.double(t(t(p)[-6] - phat[-6]) %*% obsInfoMat %*% (t(p)[-6] - phat[-6]))
      list(phat=phat,sd=sd,covar=covar,chi2 = chi2)
    })
    stopCluster(cl)
    
    res
  }


function(){
  n = 12000
  y=sim_concreteExampleFixedPartitionNoCompleteData(n)
  chi2=Reduce("c", sapply(y, function(x)x$chi2))
  phat = atrix(Reduce("c", sapply(y, function(x)x$phat)),nrow=6)
  sd = matrix(Reduce("c", sapply(y, function(x)x$sd)),nrow=6)
  covar=Reduce("+", sapply(y, function(x)x$covar))/n
  
}


function(){
  res = sim_concreteExampleFixedPartitionNoCompleteDataUnindentifiable(6000)
  apply(matrix(Reduce("c",sapply(res,function(x)x$phat)),nrow=6),1,mean) -> phat
  apply(matrix(Reduce("c",sapply(res,function(x)x$sd)),nrow=6),1,mean) -> sdmle
  q=p;print((q[2]+q[3])/(q[5]+q[6]));print((q[2]+q[5])/(q[3]+q[6]));print((q[2]+q[6])/(q[3]+q[5]))
  q=phat;print((q[2]+q[3])/(q[5]+q[6]));print((q[2]+q[5])/(q[3]+q[6]));print((q[2]+q[6])/(q[3]+q[5]))
}

function(){
  param_ident_modif1 = sim_concreteExampleFixedPartitionNoCompleteDataModif1(6000)
  apply(matrix(Reduce("c",sapply(param_ident_modif1,function(x)x$phat)),nrow=6),1,mean) -> phat
  apply(matrix(Reduce("c",sapply(param_ident_modif1,function(x)x$sd)),nrow=6),1,mean) -> sdmle
  q=p;print((q[2]+q[3])/(q[5]+q[6]));print((q[2]+q[5])/(q[3]+q[6]));print((q[2]+q[6])/(q[3]+q[5]))
  q=phat;print((q[2]+q[3])/(q[5]+q[6]));print((q[2]+q[5])/(q[3]+q[6]));print((q[2]+q[6])/(q[3]+q[5]))
}

function(){
  param_ident_modif2 = sim_concreteExampleFixedPartitionNoCompleteDataModif2(6000)
  apply(matrix(Reduce("c",sapply(param_ident_modif2,function(x)x$phat)),nrow=6),1,mean) -> phat
  apply(matrix(Reduce("c",sapply(param_ident_modif2,function(x)x$sd)),nrow=6),1,mean) -> sdmle
  q=p;print((q[2]+q[3])/(q[5]+q[6]));print((q[2]+q[5])/(q[3]+q[6]));print((q[2]+q[6])/(q[3]+q[5]))
  q=phat;print((q[2]+q[3])/(q[5]+q[6]));print((q[2]+q[5])/(q[3]+q[6]));print((q[2]+q[6])/(q[3]+q[5]))
  save("param_ident_modif2",file="E:/Dropbox/Subjects/W/Weaver Paper Reading list/Simulation Result Data/param_iden_case_modif2.RData")
}
