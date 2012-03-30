ddirtrunc <-
function(x, alpha, epsilon, log=TRUE){
	n = length(alpha)
	b = 1
	at = n*epsilon
	bt = n
	alphat = sum(alpha)
	xi = max(c( epsilon, 1-bt+b ))
	eta = min(c( b, 1-at+epsilon ))
	p = logPDFTruncatedBeta(x[n-1], alpha[n-1], alphat-alpha[n-1], xi, eta, epsilon)
	if( n > 2 ){
		for( k in (n-2):1 ){
			xi = max(c( epsilon/(1 - sum(x[(k+1):(n-1)])) , 1 - (bt - b*(n-k)) / (1-sum(x[(k+1):(n-1)])) ))
			eta = min(c( b/(1 - sum(x[(k+1):(n-1)])) , 1 - (at - epsilon*(n-k)) / (1-sum(x[(k+1):(n-1)])) ))
			tmp = x[k] / (1-sum(x[(k+1):(n-1)]))
			p = p + logPDFTruncatedBeta(tmp, alpha[k], alphat - sum(alpha[k:(n-1)]), xi, eta, epsilon)
			}
		for( k in (n-2):1 ){
			p = p - log(1 - sum(x[(k+1):(n-1)]))
			}
		}	
	if( !log ){ p = exp(p) }
	return(p)
	}

