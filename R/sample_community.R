#' Simulate sampling of a microbial community.
#'
#' `sample_community` simulates a sequencing-based sampling experiment on a
#' microbial community, with independent sequencing efforts for rDNA genes and
#' rRNA transcripts.  community. This allows for the effects of partial
#' sequencing to be evaluated.
#'
#' @return
#'
#' A tibble containing the measured rDNA and rRNA abundances drawn from the
#' community. All OTUs in the input community, including those not detected,
#' will be included. The `phantom` column indicates 'phantom OTUs' that yielded
#' rRNA transcripts but not rDNA gene sequences.
#'
#' @param community A tibble containing a microbial community. This will usually be the output of `generate_community`.
#' @param nDNA Number of rDNA gene sequences to simulate (i.e. rDNA read count)
#' @param nRNA Number of rRNA transcript sequences to simulate (i.e. rRNA read
#' count)
#'
#' @seealso generate_community
#'
#' @examples
#'
#' community <- generate_community(1000)
#' sample_community(community, nDNA = 500, nRNA = 500)
#'
#' @import tibble
#' @import dplyr
#' @export
sample_community <- function(community, nDNA = 20000, nRNA = 20000) {

  # Set up sampled community
  sampled_comm <- community %>%
    select(OTU, met_state, ribo_amp)

  # Simulate random sequencing of rDNA genes
  sampled_comm <- sample(community$OTU, nDNA, prob = community$rDNA_abund, replace = T) %>%
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
  sampled_comm <- sample(community$OTU, nRNA, prob = community$rRNA_abund, replace = T) %>%
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

  # Calculate rRNA:rDNA ratio, converting Inf (rRNA/0) and NaN (0/0) to NA
  sampled_comm <- sampled_comm %>%
    mutate(ratio = rRNA_relabund / rDNA_relabund) %>%
    mutate(ratio = case_when(
      is.na(ratio) ~ as.double(NA),
      ratio == Inf ~ as.double(NA),
      TRUE ~ ratio
    ))

  # Identify phantom OTUs (OTUs where an rRNA transcript is identified, but not
  # an rDNA gene)
  sampled_comm <- sampled_comm %>%
    mutate(phantom = rDNA_abund > 0 & rRNA_abund == 0)

  # Return
  sampled_comm
}

# Account for pesky global variables so check doesn't throw a note
globalVariables(c(
  "OTU",
  "met_state",
  "."
))
