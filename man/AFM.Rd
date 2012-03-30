\name{AFM}
\Rdversion{1.1}
\alias{AFM}
\title{
Function to fit the Admixture F-model (AFM) by MCMC
}
\description{
The function takes in genotype matrix, and returns (1) posterior of lineage loadings kappa and (2) posterior of drift parameters log(alpha).
}
\usage{
AFM(dat, nMC, burnin, thinning, eps=7, priq=1, pria=c(1,2), prik=1)
}

\arguments{
  \item{dat}{
data matrix, format as in dataset specimen
}
  \item{nMC}{
total number of MCMC iterations
}
  \item{burnin}{
no. of burnin iterations to be cast away, adjustment of proposal distributions is done during burnin
}
  \item{thinning}{
interval to sample from the MCMC, e.g. thinning=10 samples every 10th iteration
}
  \item{eps}{
truncation threshold of the Dirichlet distributions. eps=7 truncates at 1e-07, eps=8 at 1e-08, etc.
}

  \item{priq}{
specifies prior for ancestral allele frequencies, Dirichlet(priq) 
}

  \item{pria}{
specifies a prior for log(alpha): first component mean, second variance
}

  \item{prik}{
specifies a prior for kappa
} }

\details{
Convergence of the algorithm needs to be checked by the investigator (e.g. by doing postrior diagnostic tests or running a number of parallel chains.) As for the sample size, nmc=20000, burnin=10^4, thinning=10 would typically be enough, but that depends on the case. The code can be fairly slow for highly polymorphic data such as microsatellites, so try short test runs. Results may require running the code overnight. There has to be at least two subpopulations in the data. (Insert a dummy row, if you want to estimate within-population theta for one subpopulation.)

The default priors should be fairly uninformative. priq=1 implies a uniform prior, higher values indicating more symmetric frequencies. pria=c(1,2) should be flexible enough for all biologically relevant scenarios. prik scales the diagonal-dominated Dirichlet prior of kappa, higher values implying more information.
}

\value{
  \item{kapm}{Posterior of lineage loadings kappa. Dimension is (population x lineage x mc iteration).}
  \item{alpham}{Posterior of lineage-specific drift parameters log(alpha). Dimension is (lineage x mc iteration). Note that this is posterior of LOG(alpha), thus negative values are ok.}
  }

\author{
Markku Karhunen, Uni. Helsinki
}

\examples{
data(specimen)
samp <- AFM(specimen, 100, 50, 2) # few iterations!
postkap <- samp[[1]] # posterior of lineage loadings
hist(postkap[1,1,]) # histogram from population 1's loading to lineage 1
hist(postkap[1,2,]) # histogram from population 1's loading to lineage 2
}
