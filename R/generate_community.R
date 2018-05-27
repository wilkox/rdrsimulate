generate_community <- function(n = 5000) {

  # Set up the community
  comm <- tibble(OTU = 1:n)

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
    mutate(rDNA_abund = ceiling(rlnorm(n))) %>%
    mutate(rDNA_relabund = 100 * rDNA_abund / sum(rDNA_abund))

  # Set a ribosomal amplification for each OTU
  # This is based on the 'mix' model from Steven et al.
  comm <- comm %>%
    mutate(ribo_amp = case_when(
      met_state == "dead" ~ 1L,
      met_state == "dormant" ~ 100L,
      met_state == "maintenance" ~ sample(200:1000, n, replace = T),
      met_state == "growing" ~ sample(500:10000, n, replace = T)
    ))

  # Based on the rDNA abundance and ribosomal amplification, generate the rRNA
  # abundance and relative abundance, and the rRNA:rDNA ratio
  comm <- comm %>%
    mutate(rRNA_abund = rDNA_abund * ribo_amp) %>%
    mutate(rRNA_relabund = 100 * rRNA_abund / sum(rRNA_abund)) %>%
    mutate(ratio = rRNA_relabund / rDNA_relabund)

  # Return
  comm
}
