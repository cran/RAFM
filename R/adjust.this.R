adjust.this <-
function(rates, props, w, a, dirichlet=TRUE){
	adjustment = exp(a*(rates - 0.44))
	if( dirichlet ){
		output = props / adjustment
		output[which(output>10^3)] = 10^3
		output[which(output<10^(-3))] = 10^(-3)
		}
	else{ # only Gaussian distribution in this case
		output = props*adjustment
		}
	return(output)
	}

