l3 <-
function(logalpha_, prioralpha_){
	return(dnorm(logalpha_, prioralpha_[1], sqrt(prioralpha_[2]), log=TRUE))
	}

