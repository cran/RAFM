RandomDirichlet <-
function(a){
	eps = 10^(-7)
	x = rgamma(length(a), a)
	if( max(x) < 10^(-100) ){
		x = rep(0, length(a))
		x[sample(1:length(x))] = 1
		}
	x = x / sum(x)
	x[which(x<eps)] = eps
	x[which(x>1-eps)] = 1 - eps
	x = x / sum(x)
	return(x)
	}

