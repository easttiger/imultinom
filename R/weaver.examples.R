## Special input converters
Weaver.Examples.Hunter2004PL.ConvertInput <-
function(){
  #dat = read.table(file="http://sites.stat.psu.edu/~dhunter/code/btmatlab/nascar2002.txt", header=T)
  data(nascar2002_rawdata)
  dat = nascar2002_rawdata

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
	rslt$De = t(tDe);
	rslt
}

Weaver.Examples.Hankin2010Volleyball.ConvertInput <-
function(){

  data(vb_synthetic_df)
  tDe = as.matrix(vb_synthetic_df[,1:9])
  b = vb_synthetic_df$powers

  ida = tDe %*% rep(1,9) == 1L;
  a = rev(b[ida]);
  idb = !ida;

  idb[1] = FALSE; # the row has all 0
  idb[length(idb)] = FALSE; # the row has all 1
  idb[b == 0] = FALSE; # remove those b=0 cases
  b = b[idb];
  tDe = tDe[idb,];

  list(a=a,b=b,De=t(tDe))
}
