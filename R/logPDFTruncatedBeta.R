logPDFTruncatedBeta <-
function(x, alpha, beta, a, b, eps){
	if( x < a || x > b ){
		f = -Inf
		}
	else{	
		f = (alpha-1)*log(x) + (beta-1)*log(1-x) - logBeta(c(alpha,beta))
		normalizing = log(pbeta(b, alpha, beta) - pbeta(a, alpha, beta))
		if( normalizing == -Inf ){ # no place for x to move between a and b
			f = 0
			}
		else{
			f = f - normalizing
			}
		}
	return(f)
	}

