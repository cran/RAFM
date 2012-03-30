logBeta <-
function(a){
	return( sum(lgamma(a)) - lgamma(sum(a)) )
	}

