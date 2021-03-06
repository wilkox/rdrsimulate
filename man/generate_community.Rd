% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_community.R
\name{generate_community}
\alias{generate_community}
\title{Generate a simulated microbial community.}
\usage{
generate_community(n = 5000)
}
\arguments{
\item{n}{Number of OTUs in the simulated community (i.e. richness)}
}
\value{
A tibble containing the simulated community.
}
\description{
Following the approach of Steven et al. 2017, \code{generate_community} generates
a simulated microbial community, including both cell and ribosomal RNA
abundances.
}
\details{
rDNA abundances (i.e. cell counts) for each OTU are drawn from a log-normal
distribution with parameters \eqn{\mu = 0}, \eqn{\sigma = 1}.

The community is randomly divided into metabolic states of dead, dormant,
stationary or growing, with the proportion of each state drawn at random
from a uniform distribution. The metabolic state of each OTU is independent
of the OTU's rDNA abundance.

Each OTU is assigned a ribosomal amplification value, representing the
number of ribosomes per cell. This value is determined by the cell's
metabolic state:
\itemize{
\item Dead cells always have 1 ribosome per cell
\item Dormant cells always have a 100 ribosomes per cell
\item Maintenance cells have a ribosomal amplification drawn at random from the
uniform distribution U(200, 1000)
\item Growing cells have a ribosomal amplification drawn at random from the
uniform distribution U(500, 10000)
}

The rRNA:rDNA ratio is calculated for each OTU.

Note that there are some key differences to the approach described in Steven
et al. 2017:
\enumerate{
\item Only integer rDNA abundance are permitted (Steven et al. allowed
fractional abundances)
\item The ribosome amplification value for each metabolic state is drawn from a
uniform distribution (Steven et al. selected at random one of three discrete
values representing low, medium and high amplification)
}
}
\examples{

generate_community(5000)

}
\references{
Steven, B., Hesse, C., Soghigian, J., Gallegos-Graves, L. V. & Dunbar, J.
Simulated rRNA/DNA Ratios Show Potential To Misclassify Active Populations
as Dormant. Appl Env Microbiol 83, e00696-17-11 (2017).
}
\seealso{
sample_community
}
