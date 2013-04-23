AFM <-
function(dat, nMC, burnin, thinning, eps=7, priq=1, pria=c(1,2), prik=1){

	# Allele counts
	dat <- compress.data(dat)
	
	# Initial values
	epsilon <- 10^(-eps)	
	nloci <- length(dat)
	npop <- nrow(dat[[1]])
	pAnc <- dat
	z <- dat
	logalpha <- rnorm(npop, 1, sqrt(2))
	for( j in 1:nloci ){
		nal <- ncol(dat[[j]])
		pAnc[[j]] <- rdirtrunc(rep(1,nal), epsilon)
		for( i in 1:npop ){
			z[[j]][i,] <- rdirtrunc(pAnc[[j]]*exp(logalpha[i]), epsilon)
			}
		}
	kap <- matrix(0.2/(npop-1), ncol=npop, nrow=npop)
	diag(kap) <- 0.8
		
	# Initial proposals
	propAnc <- rep(10, nloci)
	propz <- matrix(100, nrow=npop, ncol=nloci)
	propkap <- rep(100, npop)
	propalpha <- rep(0.1, npop)
	iterno <- 0 # needed for adaptation only
	
	# Prior parameters
	prioralpha <- pria
	priorkap <- matrix(0.2/(npop-1), npop, npop)
	diag(priorkap) <- 0.8
	priorkap <- prik*priorkap
	priorAnc <- pAnc
	for( j in 1:nloci ){
		priorAnc[[j]] <- rep(priq, length(priorAnc[[j]]))
		}	

	# Initial likelihood
	clike1 <- rep(0, nloci)
	clike2 <- matrix(0, ncol=nloci, nrow=npop)
	clike3 <- rep(0, npop)
	clike4 <- rep(0, npop)
	clike5 <- matrix(0, ncol=nloci, nrow=npop)
	for( j in 1:nloci ){
		clike1[j] <- l1(pAnc[[j]], priorAnc[[j]], epsilon)
		for( i in 1:npop ){
			clike2[i,j] <- l2(z[[j]][i,], logalpha[i], pAnc[[j]], epsilon)
			clike5[i,j] <- l5(kap[i,], z[[j]], dat[[j]][i,])
			}
		}
	for( i in 1:npop ){
		clike3[i] <- l3(logalpha[i], prioralpha)
		clike4[i] <- l4(kap[i,], priorkap[i,], i, epsilon)
		}

	# Accept ratios
	acAnc <- rep(0, nloci)
	acz <- matrix(0, nrow=npop, ncol=nloci)
	ackap <- rep(0, npop)
	acalpha <- rep(0, npop)
	
	# Output variables
	output <- c(kap, logalpha)

	# Monte Carlo loop
	for( i in 1:nMC ){

    # Status
		if( round(i/100)==i/100 ){
			print(paste("iter", i), quote=F)
			}
		iterno <- iterno + 1
		adjust <- (i <= burnin)

    # Updating ancestral frequencies
  	for( j in 1:nloci ){
  		oldlike <- clike1[j] + sum(clike2[,j])
  		oldp <- pAnc[[j]]
  		newp <- rdirtrunc(propAnc[j]*oldp, epsilon)
  		nl1 <- l1(newp, priorAnc[[j]], epsilon)
  		nl2 <- rep(0, npop)
  		for( i in 1:npop ){
  			nl2[i] <- l2(z[[j]][i,], logalpha[i], newp, epsilon)
  			}
  		there <- ddirtrunc(newp, propAnc[j]*oldp, epsilon, log=T)
  		here <- ddirtrunc(oldp, propAnc[j]*newp, epsilon, log=T)
  		accept <- nl1 + sum(nl2) + here - oldlike - there
  		if( is.nan(accept) ){ accept <- -Inf }
  		if( log(runif(1,0,1)) < accept ){
  			pAnc[[j]] <- newp
  			clike1[j] <- nl1
  			clike2[,j] <- nl2
  			acAnc[j] <- acAnc[j] + 1
  			}
  		}

		# Updating kappa
  	for( i in 1:npop ){
  		oldkap <- kap[i,]
  		newkap <- rdirtrunc(oldkap*propkap[i], epsilon)
  		oldlike <- clike4[i] + sum(clike5[i,])
  		nl4 <- l4(newkap, priorkap[i,], i, epsilon)
  		nl5 <- rep(0, nloci)
  		for( j in 1:nloci ){
  			nl5[j] <- l5(newkap, z[[j]], dat[[j]][i,])
  			}
  		here <- ddirtrunc(oldkap, propkap[i]*newkap, epsilon)
  		there <- ddirtrunc(newkap, propkap[i]*oldkap, epsilon)
  		accept <- nl4 + sum(nl5) + here - oldlike - there
  		if( is.nan(accept) ){ accept <- -Inf }
  		if( log(runif(1,0,1)) < accept ){
  			kap[i,] <- newkap
  			clike4[i] <- nl4
  			clike5[i,] <- nl5
  			ackap[i] <- ackap[i] + 1
  			}
  		}		

    # Updating alpha
  	for( i in 1:npop ){
  		olda <- logalpha[i]
  		newa <- rnorm(1, olda, propalpha[i])
  		oldlike <- clike3[i] + sum(clike2[i,])
  		nl3 <- l3(newa, prioralpha)
  		nl2 <- rep(0, nloci)
  		for( j in 1:nloci ){
  			nl2[j] <- l2(z[[j]][i,], newa, pAnc[[j]], epsilon)
  			}
  		accept <- nl3 + sum(nl2) - oldlike # symmetric
  		if( is.nan(accept) ){
  			accept <- -Inf 
  			}
  		if( log(runif(1,0,1)) < accept ){
  			logalpha[i] <- newa
  			clike2[i,] <- nl2
  			clike3[i] <- nl3
  			acalpha[i] <- acalpha[i] + 1
  			}
  		}
  		
    # Updating frequencies in lineages
  	for( j in 1:nloci ){
  		for( i in 1:npop ){
  			oldlike <- clike2[i,j] + sum(clike5[,j])
  			oldz <- z[[j]][i,]
  			newz <- rdirtrunc(propz[i,j]*oldz, epsilon)
  			nl2 <- l2(newz, logalpha[i], pAnc[[j]], epsilon)
  			nl5 <- rep(0, npop)
  			Z <- z[[j]]
  			Z[i,] <- newz
  			for( k in 1:npop ){
  				nl5[k] <- l5(kap[k,], Z, dat[[j]][k,])
  				}
  			here <- ddirtrunc(oldz, propz[i,j]*newz, epsilon)
  			there <- ddirtrunc(newz, propz[i,j]*oldz, epsilon)
  			accept <- nl2 + sum(nl5) + here - oldlike - there
  			if( is.nan(accept) ){ accept <- -Inf }
  			if( log(runif(1,0,1)) < accept ){
  				z[[j]][i,] <- newz
  				clike2[i,j] <- nl2
  				clike5[,j] <- nl5
  				acz[i,j] <- acz[i,j] + 1
  				}
  			}
  		}
		output <- rbind(output, c(kap, logalpha)) # output variable
		
    # Adjusting proposals
  	if( adjust ){
  		w <- 1 - 0.1*exp(-i/1000)
  		a <- 0.1*exp(-i/1000)
  		iterno <- w*iterno
  		ackap <- w*ackap
  		acAnc <- w*acAnc
  		acz <- w*acz
  		acalpha <- w*acalpha
  		ratekap <- ackap / iterno
  		rateAnc <- acAnc / iterno
  		ratez <- acz / iterno
  		ratealpha <- acalpha / iterno
  		propkap <- adjust.this(ratekap, propkap, w, a)
  		propAnc <- adjust.this(rateAnc, propAnc, w, a)
  		propz <- adjust.this(ratez, propz, w, a)
  		propalpha <- adjust.this(ratealpha, propalpha, w, a, dirichlet=FALSE)
  		}	
		} # MC loop closes

	# Burnin & thinning
	output <- output[(burnin+1):nMC,]
	imax <- floor(nrow(output) / thinning)
	totake <- thinning*1:imax
	output <- output[totake,]

	# Data out
	nmc_ <- nrow(output)
	nc <- ncol(output)
	kapm <- array(NA, dim=c(npop, npop, nmc_))
	alpham <- array(NA, dim=c(npop, nmc_))
	for( i in 1:nmc_ ){
		kapm[,,i] <- matrix(output[i,1:(npop^2)], ncol=npop, nrow=npop)
		logalpha <- output[i,(npop^2+1):nc]
		alpham[,i] <- logalpha
		}
	print("posteriors written", quote=F)
	return(list(kapm,alpham))
	}

