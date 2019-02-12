#' Generate a simulated microbial community.
#'
#' Following the approach of Steven et al. 2017, `generate_community` generates
#' a simulated microbial community, including both cell and ribosomal RNA
#' abundances. 
#'
#' rDNA abundances (i.e. cell counts) for each OTU are drawn from a log-normal
#' distribution with parameters \eqn{\mu = 0}, \eqn{\sigma = 1}.
#'
#' The community is randomly divided into metabolic states of dead, dormant,
#' stationary or growing, with the proportion of each state drawn at random
#' from a uniform distribution. The metabolic state of each OTU is independent
#' of the OTU's rDNA abundance.
#'
#' Each OTU is assigned a ribosomal amplification value, representing the
#' number of ribosomes per cell. This value is determined by the cell's
#' metabolic state:
#' * Dead cells always have 1 ribosome per cell
#' * Dormant cells always have a 100 ribosomes per cell
#' * Maintenance cells have a ribosomal amplification drawn at random from the
#' uniform distribution U(200, 1000)
#' * Growing cells have a ribosomal amplification drawn at random from the
#' uniform distribution U(500, 10000)
#'
#' The rRNA:rDNA ratio is calculated for each OTU.
#'
#' Note that there are some key differences to the approach described in Steven
#' et al. 2017:
#' 1. Only integer rDNA abundance are permitted (Steven et al. allowed
#' fractional abundances)
#' 2. The ribosome amplification value for each metabolic state is drawn from a
#' uniform distribution (Steven et al. selected at random one of three discrete
#' values representing low, medium and high amplification) 
#'
#' @return
#' A tibble containing the simulated community.
#'
#' @param n Number of OTUs in the simulated community (i.e. richness)
#'
#' @seealso sample_community
#'
#' @examples
#'
#' generate_community(5000)
#'
#' @references
#'
#' Steven, B., Hesse, C., Soghigian, J., Gallegos-Graves, L. V. & Dunbar, J.
#' Simulated rRNA/DNA Ratios Show Potential To Misclassify Active Populations
#' as Dormant. Appl Env Microbiol 83, e00696-17-11 (2017).
#'
#' @import tibble
#' @import dplyr
#' @importFrom magrittr %>%
#' @export
generate_community <- function(n = 5000) {

  # Set up the community
  comm <- tibble::tibble(OTU = 1:n)

  # Randomly assign each OTU a metabolic state 
  met_distribution <- stats::runif(4)
  comm <- comm %>%
    mutate(met_state = sample(
      c("dead", "dormant", "stationary", "growing"),
      n,
      replace = T,
      prob = met_distribution
    ))

  # Randomly assign each OTU an rDNA abundance and relative abundance, based on a
  # log-normal distribution with μ = 0, σ = 1
  comm <- comm %>%
    mutate(rDNA_abund = as.integer(ceiling(stats::rlnorm(n)))) %>%
    mutate(rDNA_relabund = 100 * rDNA_abund / sum(rDNA_abund))

  # Set ribosomal amplification for each OTU
  # This is based on the 'mix' model from Steven et al.
  comm <- comm %>%
    mutate(ribo_amp = case_when(
      met_state == "dead" ~ 1L,
      met_state == "dormant" ~ 100L,
      met_state == "stationary" ~ sample(200:1000, n, replace = T),
      met_state == "growing" ~ sample(500:10000, n, replace = T)
    ))

  # Based on the rDNA abundance and ribosomal amplification, generate the rRNA
  # abundance and relative abundance, and the rRNA:rDNA ratio, converting Inf
  # (rRNA/0) and NaN (0/0) to NA
  comm <- comm %>%
    mutate(rRNA_abund = rDNA_abund * ribo_amp) %>%
    mutate(rRNA_relabund = 100 * rRNA_abund / sum(rRNA_abund)) %>%
    mutate(ratio = rRNA_relabund / rDNA_relabund) %>%
    mutate(ratio = case_when(
      is.na(ratio) ~ as.double(NA),
      ratio == Inf ~ as.double(NA),
      TRUE ~ ratio
    ))

  # Return
  comm
}

# Account for pesky global variables so check doesn't throw a note
globalVariables(c(
  "rDNA_abund",
  "ribo_amp",
  "rRNA_abund",
  "rRNA_relabund",
  "rDNA_relabund"
))
