gen.theta <-
function(kapm, alpham){

	# Eq. 12, Karhunen & Ovaskainen
	alpham = exp(alpham) # away from log-normal scale
	fstm = rep(NA, ncol(alpham))
	thetam = array(NA, dim=dim(kapm))
	npop = nrow(alpham)
	for( i in 1:length(fstm) ){
		kap = kapm[,,i]
		alpha = alpham[,i]
		theta = matrix(NA, nrow=npop, ncol=npop)
		for( j in 1:npop ){
			for( k in 1:j ){
				theta[k,j] = sum( kap[k,]*kap[j,] / (alpha + 1) )
				theta[j,k] = theta[k,j]
				}
			}
		thetam[,,i] = theta
		fstm[i] = gen.fst(theta)
		}

	return(list(thetam, fstm))
	}

