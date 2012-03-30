\name{RAFM-package}
\Rdversion{1.1}
\alias{RAFM-package}
\alias{RAFM}
\docType{package}
\title{ Admixture F-model by Metropolis-Hastings
	}

\description{ This package involves codes needed for fitting the AFM from command line, conversion of AFM to coancestry and FST, and three example datasets.

	}
\details{
\tabular{ll}{
Package: \tab RAFM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-03-30\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
For fitting the model, use program AFM. It returns (1) posterior of kappa and (2) posterior of log(alpha).
For calculating coancestry, use program gen.theta. It returns (1) posterior of coancestry matrix and (2) posterior of FST.
For data format, see e.g. data specimen.

Please note, 1: at present, other components apart from the above are undocumented. The code should be fairly transparent, and you can change it at your own risk.
Please note, 2: the convergence of the MCMC algorithm needs to be checked by the investigator, case by case. At present, the package provides no functionality for this. Secondly, the code can be very slow for polymorphic data (such as shrews and sticklebacks), so please be patient.  
}
\author{
Markku Karhunen, <markku.karhunen@helsinki.fi>
}
\examples{
data(specimen)
posterior <- AFM(specimen, 100, 50, 5) # few iterations!
theta_fst <- gen.theta(posterior[[1]], posterior[[2]])
fst <- theta_fst[[2]]
hist(fst)
fst <- gen.fst(0.8*diag(2))
fst - 0.8
}
