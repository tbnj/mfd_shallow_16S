### Setup env
#library(tidyverse)
library(dplyr)
library(stringr)

setwd("/mfd_shallow_16S/scripts/")

### Import hmm results
forward_files <- list.files("/R/hmm_forward_out", pattern = ".txt", full.names = T)
reverse_files <- list.files("/R/hmm_reverse_out", pattern = ".txt", full.names = T)

colnames <- c("target_name", "accession", "query_name", "model", "hmm_from", "hmm_to", "ali_from", 
              "ali_to", "env_from", "env_to", "sq_len", "strand", "E-value", "score", "bias", "description_of_target")

## Forward files
forward <- do.call(rbind, Map(cbind, lapply(forward_files, data.table::fread, sep = "\t", header = FALSE, fill = TRUE), 
                              file_name = basename(forward_files))) %>%
  mutate(across(V1, ~str_replace_all(., "\\s\\s*", ","))) %>%
  filter(!str_detect(V1, '#')) %>%
  tidyr::separate(V1, into = colnames, sep = ",") %>%
  mutate(across(c(hmm_from:sq_len, `E-value`:bias), ~as.numeric(.))) %>%
  #filter(hmm_from >= 27 & hmm_to <= 1391) %>%
  #select(target_name, score, file_name) %>%
  mutate(across(score, ~as.numeric(.))) %>%
  mutate(target_name = str_remove(target_name, "\\s.*"),
         source = str_remove(file_name, "_.*"))

## Remove controls
forward_reads.samples <- forward %>%
  filter(!str_detect(file_name, paste0(c("C3", "F3", "C10", "F10", "D12", "E12", "F12"), collapse = "|")))

## Reverse files
reverse <- do.call(rbind, Map(cbind, lapply(reverse_files, data.table::fread, sep = "\t", header = FALSE, fill = TRUE), 
                              file_name = basename(reverse_files))) %>%
  mutate(across(V1, ~str_replace_all(., "\\s\\s*", ","))) %>%
  filter(!str_detect(V1, '#')) %>%
  tidyr::separate(V1, into = colnames, sep = ",") %>%
  mutate(across(c(hmm_from:sq_len, `E-value`:bias), ~as.numeric(.))) %>%
  #filter(hmm_from >= 27 & hmm_to <= 1391) %>%
  #select(target_name, score, file_name) %>%
  mutate(across(score, ~as.numeric(.))) %>%
  mutate(target_name = str_remove(target_name, "\\s.*"),
         source = str_remove(file_name, "_.*"))

## Remove controls
reverse_reads.samples <- reverse %>%
  filter(!str_detect(file_name, paste0(c("C3", "F3", "C10", "F10", "D12", "E12", "F12"), collapse = "|")))

## Combine all
df.all <- rbind(forward, reverse)

## Combine without controls
all_reads.samples <- rbind(forward_reads.samples, reverse_reads.samples) %>%
  select(target_name) %>%
  distinct()

rm(forward, reverse)
gc()

## Filter for best hmm hit
df.all_filtered <- df.all %>%
  group_by(target_name) %>%
  filter(score == max(score)) %>%
  mutate(duplicate = n()) %>%
  ungroup()

rm(df.all)

## Reads found in only forward or reverse sequencing file
single <- df.all_filtered %>%
  filter(duplicate == 1)
  #select(target_name, source)

## Reads found by more than two hmms
undetermined_all_models <- df.all_filtered %>%
  filter(duplicate > 2)

## Reads found by two hmms and having equal score
undetermined_source <- df.all_filtered %>%
  filter(duplicate == 2) %>%
  select(target_name, score, source) %>%
  distinct() %>%
  group_by(target_name) %>%
  filter(2 == n())

## Filter reads found in both forward and reverse sequencing files, remove unambigious model hits and reduce to one hit only
multi <- df.all_filtered %>%
  filter(duplicate == 2,
         !target_name %in% c(undetermined_all_models, undetermined_source)) %>%
  group_by(target_name) %>%
  slice_head(n = 1)


### Filter for hits within individual regions
## ARC between 20F and SSU1000ArR
arc_reads <- single %>%
  filter(source == "arc") %>%
  rbind(multi %>% filter(source == "arc")) %>%
  filter(hmm_from >= 20 & hmm_to <= 1000) %>%
  select(target_name)

## Extract ARC read names
arc_reads.samples <- all_reads.samples %>%
  filter(target_name %in% arc_reads$target_name)

## BAC between 27F and 1391R
bac_reads <- single %>%
  filter(source == "bac") %>%
  rbind(multi %>% filter(source == "bac")) %>%
  filter(hmm_from >= 27 & hmm_to <= 1391) %>%
  select(target_name)

## Extract BAC read names
bac_reads.samples <- all_reads.samples %>%
  filter(target_name %in% bac_reads$target_name)

## EUK all hits
euk_reads <- single %>%
  filter(source == "euk") %>%
  rbind(multi %>% filter(source == "euk")) %>%
  select(target_name) %>%
  group_by(target_name)

## Extract EUK read names
euk_reads.samples <- all_reads.samples %>%
  filter(target_name %in% euk_reads$target_name)

### Sanity check for number of reads
nrow(single)+nrow(multi)

rbind(single, multi) %>% distinct() %>% nrow()

nrow(arc_reads)+nrow(bac_reads)+nrow(euk_reads)

rbind(arc_reads, bac_reads, euk_reads) %>% distinct() %>% nrow()

intersect(arc_reads, bac_reads)
intersect(arc_reads, euk_reads)
intersect(bac_reads, euk_reads)

nrow(all_reads.samples) # should be more than the next line

nrow(arc_reads.samples)+nrow(bac_reads.samples)+nrow(euk_reads.samples)

rbind(arc_reads.samples, bac_reads.samples, euk_reads.samples) %>% distinct() %>% nrow()


### Write outputs
out_path <- paste0(getwd(), "/output/")

readr::write_delim(arc_reads, paste0(out_path, format(Sys.time(), "%Y-%m-%d"), "_ARC_reads.txt"), col_names = FALSE, delim = ",")
readr::write_delim(bac_reads, paste0(out_path, format(Sys.time(), "%Y-%m-%d"), "_BAC_reads.txt"), col_names = FALSE, delim = ",")
readr::write_delim(euk_reads, paste0(out_path, format(Sys.time(), "%Y-%m-%d"), "_EUK_reads.txt"), col_names = FALSE, delim = ",")

readr::write_delim(arc_reads.samples, paste0(out_path, format(Sys.time(), "%Y-%m-%d"), "_ARC_reads_samples.txt"), col_names = FALSE, delim = ",")
readr::write_delim(bac_reads.samples, paste0(out_path, format(Sys.time(), "%Y-%m-%d"), "_BAC_reads_samples.txt"), col_names = FALSE, delim = ",")
readr::write_delim(euk_reads.samples, paste0(out_path, format(Sys.time(), "%Y-%m-%d"), "_EUK_reads_samples.txt"), col_names = FALSE, delim = ",")

rm(df.all_filtered, multi, single)
gc()
