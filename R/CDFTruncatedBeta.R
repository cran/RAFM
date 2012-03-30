CDFTruncatedBeta <-
function(x, alpha, beta, a, b){
	F = (pbeta(x, alpha, beta) - pbeta(a, alpha, beta)) / (pbeta(b, alpha, beta) - pbeta(a, alpha, beta))
	return(F)
	}

