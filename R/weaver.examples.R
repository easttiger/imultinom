
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
