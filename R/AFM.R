AFM <-
function(dat, nMC, burnin, thinning, eps=7, priq=1, pria=c(1,2), prik=1){

	# Allele counts
	dat = compress.data(dat)
	
	# Initial values
	epsilon <<- 10^(-eps)	
	nloci <<- length(dat)
	npop <<- nrow(dat[[1]])
	pAnc <<- dat
	z <<- dat
	logalpha <<- rnorm(npop, 1, sqrt(2))
	for( j in 1:nloci ){
		nal <<- ncol(dat[[j]])
		pAnc[[j]] <<- rdirtrunc(rep(1,nal), epsilon)
		for( i in 1:npop ){
			z[[j]][i,] <<- rdirtrunc(pAnc[[j]]*exp(logalpha[i]), epsilon)
			}
		}
	kap <<- matrix(0.2/(npop-1), ncol=npop, nrow=npop)
	diag(kap) <<- 0.8
		
	# Initial proposals
	propAnc <<- rep(10, nloci)
	propz <<- matrix(100, nrow=npop, ncol=nloci)
	propkap <<- rep(100, npop)
	propalpha <<- rep(0.1, npop)
	iterno <<- 0 # needed for adaptation only
	
	# Prior parameters
	prioralpha <<- pria
	priorkap <<- matrix(0.2/(npop-1), npop, npop)
	diag(priorkap) <<- 0.8
	priorkap <<- prik*priorkap
	priorAnc <<- pAnc
	for( j in 1:nloci ){
		priorAnc[[j]] <<- rep(priq, length(priorAnc[[j]]))
		}	

	# Initial likelihood
	clike1 <<- rep(0, nloci)
	clike2 <<- matrix(0, ncol=nloci, nrow=npop)
	clike3 <<- rep(0, npop)
	clike4 <<- rep(0, npop)
	clike5 <<- matrix(0, ncol=nloci, nrow=npop)
	for( j in 1:nloci ){
		clike1[j] <<- l1(pAnc[[j]], priorAnc[[j]], epsilon)
		for( i in 1:npop ){
			clike2[i,j] <<- l2(z[[j]][i,], logalpha[i], pAnc[[j]], epsilon)
			clike5[i,j] <<- l5(kap[i,], z[[j]], dat[[j]][i,])
			}
		}
	for( i in 1:npop ){
		clike3[i] <<- l3(logalpha[i], prioralpha)
		clike4[i] <<- l4(kap[i,], priorkap[i,], i, epsilon)
		}

	# Accept ratios
	acAnc <<- rep(0, nloci)
	acz <<- matrix(0, nrow=npop, ncol=nloci)
	ackap <<- rep(0, npop)
	acalpha <<- rep(0, npop)
	
	# Output variables
	output = c(kap, logalpha)

	# Monte Carlo loop
	for( i in 1:nMC ){
		if( round(i/100)==i/100 ){
			print(paste("iter", i), quote=F)
			}
		iterno <<- iterno + 1
		adjust = (i <= burnin)
		updateanc(dat, adjust)
		updatekap(dat, adjust)
		updatealpha(dat, adjust)
		updatez(dat, adjust)
		output = rbind(output, c(kap, logalpha)) # output variable
		adjust.prop(i, adjust, burnin)
		} # MC loop closes

	# Burnin & thinning
	output = output[(burnin+1):nMC,]
	imax = floor(nrow(output) / thinning)
	totake = thinning*1:imax
	output = output[totake,]

	# Data out
	nmc_ = nrow(output)
	nc = ncol(output)
	kapm = array(NA, dim=c(npop, npop, nmc_))
	alpham = array(NA, dim=c(npop, nmc_))
	for( i in 1:nmc_ ){
		kapm[,,i] = matrix(output[i,1:(npop^2)], ncol=npop, nrow=npop)
		logalpha = output[i,(npop^2+1):nc]
		alpham[,i] = logalpha
		}
	print("posteriors written", quote=F)
	return(list(kapm,alpham))
	}

