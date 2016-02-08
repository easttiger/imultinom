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
