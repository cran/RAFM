updateanc <-
function(dat, adjust){
	for( j in 1:nloci ){
		oldlike = clike1[j] + sum(clike2[,j])
		oldp = pAnc[[j]]
		newp = rdirtrunc(propAnc[j]*oldp, epsilon)
		nl1 = l1(newp, priorAnc[[j]], epsilon)
		nl2 = rep(0, npop)
		for( i in 1:npop ){
			nl2[i] = l2(z[[j]][i,], logalpha[i], newp, epsilon)
			}
		there = ddirtrunc(newp, propAnc[j]*oldp, epsilon, log=T)
		here = ddirtrunc(oldp, propAnc[j]*newp, epsilon, log=T)
		accept = nl1 + sum(nl2) + here - oldlike - there
		if( is.nan(accept) ){ accept = -Inf }
		if( log(runif(1,0,1)) < accept ){
			pAnc[[j]] <<- newp
			clike1[j] <<- nl1
			clike2[,j] <<- nl2
			acAnc[j] <<- acAnc[j] + 1
			}
		}
	}

