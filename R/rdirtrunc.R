rdirtrunc <-
function(alpha, epsilon){
	n = length(alpha)
	at = n*epsilon
	bt = n
	b = 1
	alphat = sum(alpha)
	x = rep(0,n)
	xi = max(c(epsilon, 1 - bt + b))
	eta = min(c(b, 1 - at + epsilon))
	x[n-1] = RandomTruncatedBeta(alpha[n-1], alphat-alpha[n-1], xi, eta)
	if( n > 2 ){
		for( k in (n-2):1 ){
			xi = max(c( epsilon/(1-sum(x[(k+1):(n-1)])) , 1 - (bt - b*(n-k))/(1 - sum(x[(k+1):(n-1)])) ))
			eta = min(c( b/(1-sum(x[(k+1):(n-1)])) , 1 - (at - epsilon*(n-k))/(1 - sum(x[(k+1):(n-1)])) ))
			tmp = RandomTruncatedBeta(alpha[k], alphat - sum(alpha[k:(n-1)]), xi, eta)
			x[k] = tmp*(1 - sum(x[(k+1):(n-1)]))
			}
		}
	x[n] = 1 - sum(x)
	x[which(x<epsilon)] = epsilon
	x[which(x>b)] = b
	return(x)
	}

