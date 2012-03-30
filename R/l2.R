l2 <-
function(z_, logalpha_, p_, epsilon){
	return(ddirtrunc(z_, exp(logalpha_)*p_, epsilon, log=TRUE))
	}

