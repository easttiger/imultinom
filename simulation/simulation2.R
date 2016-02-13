set.seed(1)
#configure parallel mode
options(run_parallel = T)
NCORES = as.integer(Sys.getenv("NUMBER_OF_PROCESSORS"))
require("snow")
require("imultinom")

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

      #stopCluster(cl)
      cl <- makeSOCKcluster(rep("localhost",ncluster))
      clusterExport(cl, list("p","d","R","incComplete","weaver.vanilla","weaver.bayes","Weaver"), envir = sys.frame(sys.nframe()))
      clusterExport(cl, list("randparti","mle","imultinom","g","m"), envir = .GlobalEnv)

      phat = parSapply(cl, 1:samplesize, function(x){
        lik=g(p,R,includeCompleteData = incComplete);
        if(m * 1.7^R >=d){
          weaver.vanilla(lik$a,lik$b,t(lik$De))$x
        }else{
          Weaver(lik$a,lik$b,t(lik$De))$x
        }
      })

      stopCluster(cl)
    }else{
      phat = sapply(1:samplesize, function(x){
        lik=g(p,R, includeCompleteData = incComplete);
        if(m * 1.7^R >=d){
          weaver.vanilla(lik$a,lik$b,t(lik$De))$x
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

general.partition.simulation <-
function(){
  # 4 x 4
  dev.new();
  set.seed(1)
  p=1:100;p=p/sum(p);
  par(mfrow=c(4,4));
  par(mar=c(4,5.1,1.5,1.5))
  par(mgp=c(3.5,0.6,0))
  j = 1
  t = 1:16
  s = 1:16
  for(m in c(1000,100,10,1)){
    for(R in c(16,8,4,2)){
      s[j] = m * R
      t[j] = system.time({phat=mleSim(2,p,R,2)})[["sys.self"]]
      j = j + 1
      plot(p,phat[,1],main=paste("m = ",m, ", R = ",R, sep=""), ylab=expression(hat(bolditalic(p))),xlab="");abline(0,1)
      title(xlab=expression(bolditalic(p)),line=2)
    }
  }
}
