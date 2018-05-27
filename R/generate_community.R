#' Note some refinements over Steven et al's approach:
#' 1. Only integer values for the rDNA abundance are permitted (Steven et al.
#' allowed for fractional abundances)
#' 2. For the 'mix' model of rRNA amplification, the copy number for each
#' metabolic state is drawn from a uniform distribution rather than from one of
#' three discrete values (low, medium or high amplification)
#'
#' @importFrom magrittr "%>%"
generate_community <- function(n = 5000) {

  # Set up the community
  comm <- tibble::tibble(OTU = 1:n)

  # Randomly assign each OTU a metabolic state 
  comm <- comm %>%
    mutate(met_state = sample(
      c("dead", "dormant", "maintenance", "growing"),
      n,
      replace = T
    ))

  # Randomly assign each OTU an rDNA abundance and relative abundance, based on a
  # log-normal distribution with μ = 0, σ = 1
  comm <- comm %>%
    mutate(rDNA_abund = as.integer(ceiling(rlnorm(n)))) %>%
    mutate(rDNA_relabund = 100 * rDNA_abund / sum(rDNA_abund))

  # Set ribosomal amplification for each OTU
  # This is based on the 'mix' model from Steven et al.
  comm <- comm %>%
    mutate(ribo_amp = case_when(
      met_state == "dead" ~ 1L,
      met_state == "dormant" ~ 100L,
      met_state == "maintenance" ~ sample(200:1000, n, replace = T),
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
