gen.fst <-
function(theta){
	npop = ncol(theta)
	offdiag = (sum(theta) - sum(diag(theta))) / (npop^2 - npop)
	fst = (mean(diag(theta)) - offdiag) / (1 - offdiag)
	return(fst)
	}

