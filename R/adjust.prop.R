adjust.prop <-
function(i, adjust, burnin){
	if( adjust ){
		w = 1 - 0.1*exp(-i/1000)
		a = 0.1*exp(-i/1000)
		iterno <<- w*iterno
		ackap <<- w*ackap
		acAnc <<- w*acAnc
		acz <<- w*acz
		acalpha <<- w*acalpha
		ratekap = ackap / iterno
		rateAnc = acAnc / iterno
		ratez = acz / iterno
		ratealpha = acalpha / iterno
		propkap <<- adjust.this(ratekap, propkap, w, a)
		propAnc <<- adjust.this(rateAnc, propAnc, w, a)
		propz <<- adjust.this(ratez, propz, w, a)
		propalpha <<- adjust.this(ratealpha, propalpha, w, a, dirichlet=FALSE)
		}
	}

