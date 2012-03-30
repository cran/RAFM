updatez <-
function(dat, adjust){
	for( j in 1:nloci ){
		for( i in 1:npop ){
			oldlike = clike2[i,j] + sum(clike5[,j])
			oldz = z[[j]][i,]
			newz = rdirtrunc(propz[i,j]*oldz, epsilon)
			nl2 = l2(newz, logalpha[i], pAnc[[j]], epsilon)
			nl5 = rep(0, npop)
			Z = z[[j]]
			Z[i,] = newz
			for( k in 1:npop ){
				nl5[k] = l5(kap[k,], Z, dat[[j]][k,])
				}
			here = ddirtrunc(oldz, propz[i,j]*newz, epsilon)
			there = ddirtrunc(newz, propz[i,j]*oldz, epsilon)
			accept = nl2 + sum(nl5) + here - oldlike - there
			if( is.nan(accept) ){ accept = -Inf }
			if( log(runif(1,0,1)) < accept ){
				z[[j]][i,] <<- newz
				clike2[i,j] <<- nl2
				clike5[,j] <<- nl5
				acz[i,j] <<- acz[i,j] + 1
				}
			}
		}
	}

