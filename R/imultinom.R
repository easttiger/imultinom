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
