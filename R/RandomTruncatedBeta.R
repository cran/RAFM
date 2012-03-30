RandomTruncatedBeta <-
function(alpha, beta, a, b){
	quantil = pbeta(a, alpha, beta) + runif(1,0,1)*(pbeta(b, alpha, beta) - pbeta(a, alpha, beta))
	x = qbeta(quantil, alpha, beta)
	return(x)
	}

