logPDFDirichlet <-
function(x, a){
	return( sum((a-1)*log(x)) - logBeta(a) )
	}

