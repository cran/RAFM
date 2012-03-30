l5 <-
function(kap_, z_, dat){
	K = matrix(kap_, nrow=nrow(z_), ncol=ncol(z_))
	p = apply(K*z_, 2, sum)
	f = dmultinom(dat, size=sum(dat), prob=p, log=TRUE)
 	return(f)
	}

