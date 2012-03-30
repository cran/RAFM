updatekap <-
function(dat, adjust){
	for( i in 1:npop ){
		oldkap = kap[i,]
		newkap = rdirtrunc(oldkap*propkap[i], epsilon)
		oldlike = clike4[i] + sum(clike5[i,])
		nl4 = l4(newkap, priorkap[i,], i, epsilon)
		nl5 = rep(0, nloci)
		for( j in 1:nloci ){
			nl5[j] = l5(newkap, z[[j]], dat[[j]][i,])
			}
		here = ddirtrunc(oldkap, propkap[i]*newkap, epsilon)
		there = ddirtrunc(newkap, propkap[i]*oldkap, epsilon)
		accept = nl4 + sum(nl5) + here - oldlike - there
		if( is.nan(accept) ){ accept = -Inf }
		if( log(runif(1,0,1)) < accept ){
			kap[i,] <<- newkap
			clike4[i] <<- nl4
			clike5[i,] <<- nl5
			ackap[i] <<- ackap[i] + 1
			}
		}
	}

