#setwd('C:/temp/Real Data Examples/')
# Please manually set the current directory which should contain
#   - demo_weaver.R (this script)
#   - weaver.R
#   - 3 data folders Hankin2010_Volleyball/, Hunter2004_nascar2002_BT/, and Hunter2004_nascar2002_PL/

GLOBAL_SWITCH_RECORD_ITERATION = FALSE # try turn this to TRUE -> to see iteration details and plots

cat('Sourcing weaver.R from current directory\n')
source(file="weaver.R");

cat("\n=================\n")
subject = "Basic Weaver on Ford (1957) adjacency matrix format"
cat(c(subject,"\n"))
Ford1957_adjmat = matrix( 
	c( 0,15,15, 0,
	  11, 0,10,20,
	  11,10, 0,20,
	   0, 1, 1, 0),
	nrow = 4,
	byrow=TRUE);
Ford1957 = Weaver.input.matrix(Ford1957_adjmat);
Ford1957_out = WeaverBas(Ford1957$a, Ford1957$b, Ford1957$tDe,iteration=GLOBAL_SWITCH_RECORD_ITERATION)
if(GLOBAL_SWITCH_RECORD_ITERATION){
  with(Ford1957_out, print(list(x=x,iterCount=iter$count,lnLik=iter$lnLik[iter$count])))
  plot(Ford1957_out$iter$lnLik, type="b",pch=20,main=subject)
}else{
  print(Ford1957_out)
}

cat("\n=================\n")
subject = "Basic Weaver on Hunter (2004) NASCAR2002 Plackett-Luce Modelling"
cat(c(subject,"\n"))
a = read.table(file="Hunter2004_nascar2002_PL/a.txt",header=F)[[1]]
b = read.table(file="Hunter2004_nascar2002_PL/b.txt",header=F)[[1]]
tDe = as.matrix(read.table(file="Hunter2004_nascar2002_PL/tDe.txt",header=F))
nascar2002PL = list(a=a,b=b,tDe=tDe);
# alternatively, use following line to read data
#   nascar2002PL = Weaver.Examples.Hunter2004PL.ConvertInput
nascar2002PL_out = WeaverBas(nascar2002PL$a,nascar2002PL$b,nascar2002PL$tDe,iteration=GLOBAL_SWITCH_RECORD_ITERATION)
if(GLOBAL_SWITCH_RECORD_ITERATION){
  with(nascar2002PL_out, print(list(x=x,iterCount=iter$count,lnLik=iter$lnLik[iter$count])))
  x11();
  plot(nascar2002PL_out$iter$lnLik, type="b",pch=20,main=subject)
}else{
  print(nascar2002PL_out)
}


cat("\n=================\n")
subject = "Basic Weaver on Hunter (2004) NASCAR2002 Bradley-Terry Modelling"
cat(c(subject,"\n"))
a = read.table(file="Hunter2004_nascar2002_BT/a.txt",header=F)[[1]]
b = read.table(file="Hunter2004_nascar2002_BT/b.txt",header=F)[[1]]
tDe = as.matrix(read.table(file="Hunter2004_nascar2002_BT/tDe.txt",header=F))
nascar2002BT = list(a=a,b=b,tDe=tDe);
nascar2002BT_out = WeaverBas(nascar2002BT$a,nascar2002BT$b,nascar2002BT$tDe,iteration=GLOBAL_SWITCH_RECORD_ITERATION)
if(GLOBAL_SWITCH_RECORD_ITERATION){
  with(nascar2002BT_out, print(list(x=x,iterCount=iter$count,lnLik=iter$lnLik[iter$count])))
  x11();
  plot(nascar2002BT_out$iter$lnLik, type="b",pch=20,main=subject)
}else{
  print(nascar2002BT_out)
}


cat("\n=================\n")
subject = "Greedy Weaver on Hankin (2010) Volleyball dataset";
cat(c(subject,"\n"))
a = read.table(file="Hankin2010_Volleyball/a.txt",header=F)[[1]]
b = read.table(file="Hankin2010_Volleyball/b.txt",header=F)[[1]]
tDe = as.matrix(read.table(file="Hankin2010_Volleyball/tDe.txt",header=F))
Hankin2010 = list(a=a,b=b,tDe=tDe);
# alternatively, use following line to read data
#   Hankin2010 = Weaver.Examples.Hankin2010Volleyball.ConvertInput()
Hankin2010_out = WeaverGre(Hankin2010$a,Hankin2010$b,Hankin2010$tDe,iteration=GLOBAL_SWITCH_RECORD_ITERATION)
if(GLOBAL_SWITCH_RECORD_ITERATION){
  with(Hankin2010_out, print(list(x=x,iterCount=iter$count,lnLik=iter$lnLik[iter$count])))
  x11();
  plot(Hankin2010_out$iter$lnLik, type="b",pch=20,main=subject)
}else{
  print(Hankin2010_out)
}


cat("\n=================\n")
subject = "Bayeisan Weaver on Hankin (2010) Volleyball dataset"
cat(c(subject,"\n"))
Hankin2010_out2 = WeaverBayes(Hankin2010$a,Hankin2010$b,Hankin2010$tDe,iteration=GLOBAL_SWITCH_RECORD_ITERATION)
if(GLOBAL_SWITCH_RECORD_ITERATION){
  with(Hankin2010_out2, print(list(x=x,iterCount=iter$count,lnLik=iter$lnLik[length(iter$lnLik)])))
  x11();
  plot(Hankin2010_out2$iter$lnLik,type="l",main=subject)
}else{
  print(Hankin2010_out2)
}
