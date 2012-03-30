compress.data <-
function(rawDat){
	pops = rawDat[,1]
	DNA = rawDat[,2:ncol(rawDat)]
	dna = matrix(NA, ncol=ncol(DNA)/2, nrow=2*nrow(DNA))
	for( i in 1:(ncol(DNA)/2) ){
		dna[,i] = c(DNA[,2*i], DNA[,-1+2*i])
		}
	pops = rep(pops, 2)
	uniqs = unique(pops)
	nloci = ncol(dna)
	dat = list(1:nloci)
	for( j in 1:nloci ){
		dnaj = dna[,j]
		nas = which(is.na(dnaj))
		if( length(nas) > 0 ){
			dnaj = dnaj[-nas]
			popsj = pops[-nas]
			}
		else{ popsj = pops }
		genotypes = c(unique(dnaj), "unobs")
		datj = matrix(0, nrow=length(uniqs), ncol=length(genotypes))
		for( p in 1:length(uniqs) ){
			thispop = dnaj[which(popsj==uniqs[p])]
			datj[p, ] = sapply(genotypes, function(x) length(which(thispop==x)))
			}
		dat[[j]] = datj
		}
	return(dat)
	}

