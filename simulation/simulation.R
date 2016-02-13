set.seed(1)
#configure parallel mode
options(run_parallel = T)
NCORES = as.integer(Sys.getenv("NUMBER_OF_PROCESSORS"))
require("snow")
require("imultinom")
p = t(matrix(c(
  0.1654,0.2024,0.1444,0.1532,0.2301,0.1046
),nrow = 3))
subsize = c(120,40,40,100)
De = matrix(c(
  1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1
),nrow = 6)

sim_concreteExampleFixedPartition <-
  function(samplesize = 1) {
    # this block simulates from 1.2 A concrete example
    # fix: partition and sample sizes
    res = NULL
    im = imultinom(a = double(nrow(De)),b = double(ncol(De)),De = De)

    cl <- makeSOCKcluster(rep("localhost",NCORES))
    clusterExport(cl, list("p","subsize","De","im","mle"),envir = environment())
    res = parLapply(cl, 1:samplesize, function(x) {
      y = list()
      y[[1]] = as.vector(table(sample(
        x = c(1,4,2,5,3,6),prob = p,size = subsize[1],replace = T
      )))
      y[[2]] = as.vector(table(sample(
        x = c("135","246"),prob = apply(p,1,sum),size = subsize[2],replace = T
      )))
      y[[3]] = as.vector(table(sample(
        x = c("14","25","36"),prob = apply(p,2,sum),size = subsize[3],replace =
          T
      )))
      y[[3]][1] = y[[3]][1] - 100
      y[[4]] = as.vector(table(sample(
        x = c("1","4"),prob = p[,1],size = subsize[4],replace = T
      )))
      a = y[[1]]
      a[1] = a[1] + y[[4]][1]
      a[4] = a[4] + y[[4]][2]
      b = c(y[[2]], y[[3]])
      im@a = a;
      im@b = b;
      mle = mle(im)
      obsInfoMat = solve(mle$covmle[1:5,1:5])
      phat = mle$x
      covphat = mle$covmle
      list(
        phat = phat,
        covphat = covphat,
        sdphat = sqrt(as.vector(diag(covphat))),
        chi2 = as.double(t(t(p)[-6] - phat[-6]) %*% obsInfoMat %*% (t(p)[-6] - phat[-6]))
      )
    })
    stopCluster(cl)
    save(
      res,file = paste(
        "c:/temp/sim_concreteExampleFixedPartition(",samplesize,")",format(Sys.time(),"_%Y-%m-%d_%H.%M.%S.rda"),sep =
          ""
      )
    )
    res
  }

coverage.probability <-
  function(n) {
    set.seed(1)
    #best n=12000
    proba = seq(1.0 / n,1.0,1.0 / n)
    y = sim_concreteExampleFixedPartition(n)
    phat = matrix(nrow = 6,Reduce("c", sapply(y, function(x)
      (x$phat))))
    write.csv(t(phat),format(Sys.time(),"c:/temp/phat_%Y-%m-%d_%H.%M.%S.csv"))
    sdphat = matrix(nrow = 6,Reduce("c", sapply(y, function(x)
      (x$sdphat))))
    z = (phat - as.vector(t(p))) / sdphat
    chi2 = Reduce("c", sapply(y, function(x)
      x$chi2))
    #layout(matrix(c(1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13),nrow=6,ncol=3,byrow=T))
    par(mfrow = c(1,1))
    plot(
      proba,sapply(proba, function(x) {
        sum(chi2 < qchisq(p = x,df = 5)) / n
      }),main = "Joint Coverage Probability\n at all conficence levels", xlab =
        "Cumulative Probabilities of Chisq(5) (The Confidence Level)", ylab = "Empirical Cumulative Probabilities (1 - pvalue)", type = 'p',cex = 0.2
    )
    abline(a=0,b=1,col=2)
    par(mfrow = c(4,3))
    par(cex.axis=0.8)
    par(mar=c(2.5,2.1,1.3,0.8))
    par(mgp=c(1.4,0.5,0))
    for (i in 1:6) {
      plot(
        proba,sapply(proba, function(x) {
          sum(z[i,] < qnorm(p = x)) / n
        }), main = bquote(paste(
          "Marginal Coverage Pr. of ", p[.(i)], " all C.L.", sep = ""
        )),cex.main = 0.9, xlab = paste(
          "Cum. Pr. N(", format(100 * t(p)[i],digits = 4), "%, (", format(100 * mean(sdphat[i,]),digits = 3), "%)^2)", sep = ""
        ),ylab = "Emp.Cum.Pr.",cex = 0.1
      )
      abline(a=0,b=1,col=2)
    }

    for (i in 1:6) {
      x = phat[i,]
      hist(
        x,prob = T,nclass = 50, xlab=NULL,ylab = NULL,xaxt="n",main = bquote(paste("Marginal Emp. Dist. of ", hat(p)[.(i)])),cex.main =
          0.9
      ) #, bquote(hat(p)[.(i)])
      sd = mean(sdphat[i,])
      axis(1,at = seq(round(t(p)[i]-3.5*sd,2),round(t(p)[i]+3.5*sd,2),length.out=5
      ))
      lines(seq(min(x),max(x),0.001),dnorm(
        seq(min(x),max(x),0.001),mean = t(p)[i],sd = mean(sdphat[i,])
      ),col = 2)
    }
  }

get.mean.covphat <-
function(folder="E:/Dropbox/GitHub/imultinom/simulation/results/"){
  load(paste(folder,"sim_concreteExampleFixedPartition(60000)_2016-02-10_02.35.33.rda", sep=""))
  s = Reduce("+", lapply(res, function(x){x$covphat}))
  se = Reduce("+", lapply(res, function(x){x$sdphat}))
  n = length(res)
  rm(res)
  load(paste(folder,"sim_concreteExampleFixedPartition(60000)_2016-02-10_02.54.32.rda", sep=""))
  s = Reduce("+", lapply(res, function(x){x$covphat}), s)
  se = Reduce("+", lapply(res, function(x){x$sdphat}), se)
  n = n + length(res)
  rm(res)
  list(mean.cov = s / n, mean.se = se / n)

}
