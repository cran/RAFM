\name{gen.fst}
\Rdversion{1.1}
\alias{gen.fst}

\title{
Calculation of FST out of a coancestry matrix.
}

\description{
A very simple function of the theta matrix.
}

\usage{
gen.fst(theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
	population-to-population coancestry matrix
}
}

\details{
Weir-Cockerham definition of FST, known to be consistent with Wright and coalescent theory.
}

\value{
Scalar number in [0,1].
}

\references{
Cockerham, C. C., and B. S. Weir, 1987 Correlations, Descent Measures - Drift with Migration and Mutation. PNAS 84: 8512-8514.
}

\author{
Markku Karhunen, Uni. Helsinki
}

\examples{
theta <- 0.15*diag(3)
gen.fst(theta)
}

