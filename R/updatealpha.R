updatealpha <-
function(dat, adjust){
	for( i in 1:npop ){
		olda = logalpha[i]
		newa = rnorm(1, olda, propalpha[i])
		oldlike = clike3[i] + sum(clike2[i,])
		nl3 = l3(newa, prioralpha)
		nl2 = rep(0, nloci)
		for( j in 1:nloci ){
			nl2[j] = l2(z[[j]][i,], newa, pAnc[[j]], epsilon)
			}
		accept = nl3 + sum(nl2) - oldlike # symmetric
		if( is.nan(accept) ){
			accept = -Inf 
			}
		if( log(runif(1,0,1)) < accept ){
			logalpha[i] <<- newa
			clike2[i,] <<- nl2
			clike3[i] <<- nl3
			acalpha[i] <<- acalpha[i] + 1
			}
		}
	}

