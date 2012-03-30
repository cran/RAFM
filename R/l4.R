l4 <-
function(kap_, priorkap_, i_, epsilon){
	if( max(kap_)==kap_[i_] ){
		f = ddirtrunc(kap_, priorkap_, epsilon, log=TRUE)
		}
	else{ f = -Inf }
	return(f)	
	}

