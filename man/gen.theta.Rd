\name{gen.theta}
\Rdversion{1.1}
\alias{gen.theta}

\title{
Conversion of AFM parameters to coancestry matrix and FST.
}
\description{
The AFM parameters involve the lineage loadings (kappa) and the drift parameters (alpha). Theta is the coancestry matrix, and FST is Wright's fixation index calculated out of theta. This function does the deterministic conversion of posteriors: (kappa, alpha) -> (theta, FST).
}
\usage{
gen.theta(kapm, alpham)
}

\arguments{
  \item{kapm}{
	posterior of the kappa matrix, given by function AFM. Dimension must be (population x lineage x mc iteration).
	}
  \item{alpham}{
	posterior of the log(alpha) vector, given by the function AFM. Dimension must be (lineage x mc iteration).
}
}

\details{
Note that input 2 is the posterior of LOG(alpha), i.e. it goes on a scale from minus to plus infinity.
}

\value{
\item{ thetam }{ posterior of coancestry matrix, dimension (population x population x mc iteration) }
\item{ fstm }{ posterior of FST }
}

\author{
Markku Karhunen, Uni. Helsinki
}

\examples{
data(specimen)
post_AFM <- AFM(specimen, 10, 4, 2) # a silly posterior with too few iterations
kap <- post_AFM[[1]]
logalpha <- post_AFM[[2]]
post <- gen.theta(kap, logalpha) # conversion
n <- sample(1:3, 1)
post[[1]][,,n] # one sample from theta's posterior, should be a symmetric matrix
}

