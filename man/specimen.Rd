\name{specimen}
\Rdversion{1.1}
\alias{specimen}
\docType{data}
\title{
A small simulated data.
}
\description{
Format your own data according to this sample dataset for use in AFM. 
}
\usage{data(specimen)}
\format{
  A data frame with 100 observations on the following 21 variables.
  \describe{
    \item{\code{subpop}}{subpopulation code}
    \item{\code{loc1.1}}{allelic copy 1, locus 1}
    \item{\code{loc1.2}}{allelic copy 2, locus 1}
    \item{\code{loc2.1}}{allelic copy 1, locus 2}
    \item{\code{loc2.2}}{allelic copy 2, locus 2}
    \item{\code{loc3.1}}{allelic copy 1, locus 3}
    \item{\code{loc3.2}}{allelic copy 2, locus 3}
    \item{\code{loc4.1}}{allelic copy 1, locus 4}
    \item{\code{loc4.2}}{allelic copy 2, locus 4}
    \item{\code{loc5.1}}{allelic copy 1, locus 5}
    \item{\code{loc5.2}}{allelic copy 2, locus 5}
    \item{\code{loc6.1}}{allelic copy 1, locus 6}
    \item{\code{loc6.2}}{allelic copy 2, locus 6}
    \item{\code{loc7.1}}{allelic copy 1, locus 7}
    \item{\code{loc7.2}}{allelic copy 2, locus 7}
    \item{\code{loc8.1}}{allelic copy 1, locus 8}
    \item{\code{loc8.2}}{allelic copy 2, locus 8}
    \item{\code{loc9.1}}{allelic copy 1, locus 9}
    \item{\code{loc9.2}}{allelic copy 2, locus 9}
    \item{\code{loc10.1}}{allelic copy 1, locus 10}
    \item{\code{loc10.2}}{allelic copy 2, locus 10}
  }
}
\details{
Each row represents one individual. Column 1 contains the subpopulation code (can be any number or string). Columns 2+ contain the genotypic information; genotypes can be coded by any number or string. Each number/string is one haploid genotype, while suites of two columns represent diploid genotypes. Missing data must be coded by NA. (Missing the other half of a diploid genotype is ok.) Subpopulation codes must NOT be missing, and there has to be at least two subpopulations. (Insert a dummy row, if you want to estimate within-population theta for one subpopulation.) 
}

\examples{
data(specimen)
head(specimen) # how it looks
write.table(specimen, row.names=FALSE, file="specimen.txt") # for reference
}
\keyword{datasets}
