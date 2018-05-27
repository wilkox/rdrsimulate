# This will become the sample_community function
sample_community <- function(community, nDNA = 20000, nRNA = 20000) {

  # Set up sampled community
  sampled_comm <- comm %>%
    select(OTU, met_state, ribo_amp)

  # Simulate random sequencing of rDNA genes
  sampled_comm <- sample(comm$OTU, nDNA, prob = comm$rDNA_abund, replace = T) %>%
    tibble(OTU = .) %>%
    group_by(OTU) %>%
    summarise(rDNA_abund = n()) %>%
    ungroup() %>%
    right_join(sampled_comm, by = "OTU") %>%
    mutate(rDNA_abund = case_when(
      is.na(rDNA_abund) ~ 0L,
      TRUE ~ rDNA_abund
    )) %>%
    mutate(rDNA_relabund = 100 * rDNA_abund / sum(rDNA_abund))

  # Simulate random sequencing of rRNA transcripts
  sampled_comm <- sample(comm$OTU, nRNA, prob = comm$rRNA_abund, replace = T) %>%
    tibble(OTU = .) %>%
    group_by(OTU) %>%
    summarise(rRNA_abund = n()) %>%
    ungroup() %>%
    right_join(sampled_comm, by = "OTU") %>%
    mutate(rRNA_abund = case_when(
      is.na(rRNA_abund) ~ 0L,
      TRUE ~ rRNA_abund
    )) %>%
    mutate(rRNA_relabund = 100 * rRNA_abund / sum(rRNA_abund))

  # Reorder columns
  sampled_comm <- select(
    sampled_comm,
    OTU,
    met_state,
    rDNA_abund,
    rDNA_relabund,
    ribo_amp,
    rRNA_abund,
    rRNA_relabund
  )

  # Calculate rRNA:rDNA ratio
  sampled_comm <- sampled_comm %>%
    mutate(ratio = rDNA_relabund / rRNA_relabund)

  # Identify phantom OTUs (OTUs where an rRNA transcript is identified, but not
  # an rDNA gene)
  sampled_comm <- sampled_comm %>%
    mutate(phantom = rDNA_abund > 0 & rRNA_abund == 0)

  # Return
  sampled_comm
}
