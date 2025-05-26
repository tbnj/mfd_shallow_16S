### Setup env
#library(tidyverse)
library(dplyr)
library(stringr)

setwd("/mfd_shallow_16S")

### Import data

## Load combined sintax file
data <- data.table::fread("sintax_combined_out/DATE_arcbac_sintax_combined.csv", sep = ",", header = TRUE) %>%
  rename(read_name = Read_name)

### Import SINTAX data

## List files
forward_files <- list.files("sintax_forward_out", pattern = "arc_bac_", full.names = F)
reverse_files <- list.files("sintax_reverse_out", pattern = "arc_bac_", full.names = F)

## Import forward
df.forward <- data.frame(seq_id = forward_files) %>%
  mutate(across(seq_id, ~str_remove(., "arc_bac_")),
         across(seq_id, ~str_remove(., "_forward_MFG_ssu_database_NR987_trunc.sintax")))

## Import reverse
df.reverse <- data.frame(seq_id = reverse_files) %>%
  mutate(across(seq_id, ~str_remove(., "arc_bac_")),
         across(seq_id, ~str_remove(., "_reverse_MFG_ssu_database_NR987_trunc.sintax")))

#intersect <- intersect(tmp1, tmp2)

df.forward.filt <- df.forward %>%
  filter(seq_id %in% intersect$seq_id) %>%
  mutate(start = "arc_bac_",
         path = "sintax_forward_out/",
         end = "_forward_MFG_ssu_database_NR987_trunc.sintax") %>%
  mutate(id = str_c(path, start, seq_id, end))

df.reverse.filt <- df.reverse %>%
  filter(seq_id %in% intersect$seq_id) %>%
  mutate(start = "arc_bac_",
         path = "sintax_reverse_out/",
         end = "_reverse_MFG_ssu_database_NR987_trunc.sintax") %>%
  mutate(id = str_c(path, start, seq_id, end))

## Bind together and pull id
files <- rbind(df.forward.filt, df.reverse.filt) %>%
  arrange(seq_id) %>%
  pull(id)

## Split data 
files.list <- split(files, ceiling(seq_along(files)/100))


### Import metadata

## Sample ID to library ID
meta.samples.collapsed <- data.table::fread("metadata/2025-04-24_corrected_combined_metadata_cleaned.csv", sep = ",", header = TRUE) %>%
  select(fieldsample_barcode, library_id)

## Seq ID to library ID
meta.samples <- data.table::fread("metadata/2024-02-20_all_samples_metadata.csv", sep = ",", header = TRUE) %>%
  select(seq_id, library_id) %>%
  filter(!is.na(seq_id))

## Controls
meta.controls <- data.table::fread("metadata/2024-02-20_mfd_controls_metadata.csv", sep = ",", header = TRUE) %>%
  select(seq_id, library_id) %>%
  filter(!is.na(seq_id))

## Mapping of samples
mapping.samples <- data.table::fread("metadata/2025-04-24_corrected_combined_metadata_cleaned.csv", sep = ",", header = TRUE) %>%
  select(fieldsample_barcode, seq_id, library_id) %>%
  filter(!is.na(seq_id)) %>%
  select(fieldsample_barcode, library_id)

## Combine metadata
metadata <- meta.samples %>%
  rbind(meta.controls)


### Generate key-pairs
create.key <- function(list, metadata) {
  list <- list

  colnames <- c("read_name", "tax_string", "V4", "taxonomy", "seq_id")

  for (i in 1:length(list)) {
    key <- setNames(do.call(rbind, Map(cbind, lapply(list[[i]], data.table::fread, sep = "\t",
                                                     header = FALSE, fill = TRUE), V5 = basename(list[[i]]))),
                    colnames) %>%
      select(read_name, seq_id) %>%
      mutate(across(read_name, ~str_remove(., "\\s.*")),
             across(seq_id, ~str_remove(., "arc_bac_")),
             across(seq_id, ~str_remove(., "_forward_MFD_ssu_database_NR987_trunc.sintax")),
             across(seq_id, ~str_remove(., "_reverse_MFD_ssu_database_NR987_trunc.sintax")),
             library_id = str_remove(seq_id, "_[^_]+$")) %>%
      distinct() %>%
      filter(read_name %in% data$read_name)

    print(i)

    list[[i]] <- key
  }

  df <- bind_rows(list)

  return(df)
}

## Combine key-pairs
key <- create.key(list = files.list, metadata = metadata)

## Write to disk
data.table::fwrite(key, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_arcbac_MFD_read_key.csv"))

## Import key-pairs
# key <- data.table::fread("output/2023-10-12_arcbac_MFD_read_key.csv", sep = ",", header = TRUE)


### Create a temporary combined table

## Change data to long format
df.long  <- key %>%
  left_join(data, by = "read_name") %>%
  select(-seq_id, library_id, everything()) %>%
  select(-read_name) %>%
  tidyr::unite(col = "taxonomy", 2:8, sep = ",") %>%
  select(library_id, taxonomy) %>%
  group_by(library_id, taxonomy) %>%
  summarise(Count = n()) %>%
  ungroup()

## Split into lists
libs.list <- split(df.long, df.long$library_id)

#tmp.list <- libs.list[c(1:200)]

split.list <- split(libs.list, ceiling(seq_along(libs.list)/100))

## Define function for creating temporary tables
create.tmp.phylotable <- function(list) {
  list <- list

  list.tmp <- list()

  for (i in 1:length(list)) {
    tmp <- list[[i]] %>%
      group_by(taxonomy) %>%
      #tidyr::pivot_wider(., names_from = "library_id", values_from = "Count") %>%
      tidyr::spread(., key = "library_id", value = "Count") %>%
      ungroup()

    print(i)

    list.tmp[[i]] <- tmp
  }

  list.df <- purrr::reduce(list.tmp, full_join, by = 'taxonomy') %>%
    select(-taxonomy, taxonomy) %>%

  return(list.df)
}


## Run function
phylotable.tmp <- lapply(split.list, create.tmp.phylotable) %>%
  purrr::reduce(., full_join, by = 'taxonomy') %>%
  select(-taxonomy, taxonomy) %>%
  #mutate(across(where(is.numeric), ~tidyr::replace_na(., 0))) %>%
  mutate(OTU = paste0("OTU_", 1:nrow(.))) %>%
  select(OTU, everything())

## Make "OTU" key
otukey <- phylotable.tmp %>%
  select(OTU, taxonomy)

## Select library IDs of samples
samples <- phylotable.tmp %>%
  colnames() %>%
  str_subset(pattern = "LIB")


samples_first <- samples[c(1:(length(samples)/2))]

samples_second <- samples[c((length(samples)/2+1):length(samples))]

intersect(samples_first, samples_second)

## Divide temporary table in two due to two many rows in long format
phylotable.tmp_first <- phylotable.tmp %>%
  select(any_of(samples_first), taxonomy)

phylotable.tmp_second <- phylotable.tmp %>%
  select(any_of(samples_second), taxonomy)

## Change temporary tables to long format

### Controls
phylotable.long.controls_first <- phylotable.tmp_first %>%
  select(any_of(meta.controls$library_id), taxonomy) %>%
  group_by(taxonomy) %>%
  tidyr::pivot_longer(cols = starts_with("LIB"), names_to = "library_id", values_to = "Count") %>%
  filter(!is.na(Count))

phylotable.long.controls_second <- phylotable.tmp_second %>%
  select(any_of(meta.controls$library_id), taxonomy) %>%
  group_by(taxonomy) %>%
  tidyr::pivot_longer(cols = starts_with("LIB"), names_to = "library_id", values_to = "Count") %>%
  filter(!is.na(Count))

### Samples
phylotable.long.samples_first <- phylotable.tmp_first %>%
  select(any_of(meta.samples$library_id), taxonomy) %>%
  group_by(taxonomy) %>%
  tidyr::pivot_longer(cols = starts_with("LIB"), names_to = "library_id", values_to = "Count") %>%
  filter(!is.na(Count)) %>%
  left_join(mapping.samples, relationship = "many-to-many") %>%
  select(-library_id) %>%
  group_by(fieldsample_barcode, taxonomy) %>%
  mutate(Count_sum = sum(Count), .keep = "unused") %>%
  ungroup() %>%
  rename(Count = Count_sum) %>%
  distinct()

phylotable.long.samples_second <- phylotable.tmp_second %>%
  select(any_of(meta.samples$library_id), taxonomy) %>%
  group_by(taxonomy) %>%
  tidyr::pivot_longer(cols = starts_with("LIB"), names_to = "library_id", values_to = "Count") %>%
  filter(!is.na(Count)) %>%
  left_join(mapping.samples, relationship = "many-to-many") %>%
  select(-library_id) %>%
  group_by(fieldsample_barcode, taxonomy) %>%
  mutate(Count_sum = sum(Count), .keep = "unused") %>%
  ungroup() %>%
  rename(Count = Count_sum) %>%
  distinct()

## Create specific lists
phylotable.controls.list <- c(split(phylotable.long.controls_first, phylotable.long.controls_first$library_id),
                              split(phylotable.long.controls_second, phylotable.long.controls_second$library_id))


phylotable.samples.list <- c(split(phylotable.long.samples_first, phylotable.long.samples_first$fieldsample_barcode),
                             split(phylotable.long.samples_second, phylotable.long.samples_second$fieldsample_barcode))

## Split
split.list.controls <- split(phylotable.controls.list, ceiling(seq_along(phylotable.controls.list)/100))
split.list.samples <- split(phylotable.samples.list, ceiling(seq_along(phylotable.samples.list)/100))


## Define new function for creating final tables
create.phylotable <- function(list, string) {
  list <- list

  list.tmp <- list()

  for (i in 1:length(list)) {
    phylo <- list[[i]] %>%
      group_by(taxonomy) %>%
      tidyr::pivot_wider(., names_from = string, values_from = "Count")

    print(i)

    list.tmp[[i]] <- phylo
  }

  list.df <- purrr::reduce(list.tmp, full_join, by = 'taxonomy') %>%
    select(-taxonomy, taxonomy)

  return(list.df)
}

### Create final phylotables
## Controls
phylotable.controls <- lapply(split.list.controls, create.phylotable, "library_id") %>%
  purrr::reduce(., full_join, by = 'taxonomy') %>%
  select(-taxonomy, taxonomy) %>%
  left_join(otukey, by = 'taxonomy') %>%
  select(OTU, starts_with("LIB"), taxonomy) %>%
  arrange(., str_rank(OTU, numeric = TRUE)) %>%
  #mutate(across(where(is.numeric), ~tidyr::replace_na(., 0))) %>% # Do this with UNIX sed command
  tidyr::separate_wider_delim(taxonomy, delim = ",",
                              names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                              too_few = "align_start")

## Samples
phylotable.samples.tmp <- lapply(split.list.samples, create.phylotable, "fieldsample_barcode") %>%
  purrr::reduce(., full_join, by = 'taxonomy') %>%
  select(-taxonomy, taxonomy)

## Fix duplex columns
colnames.dup <- phylotable.samples.tmp %>%
  select(starts_with("MFD"), taxonomy) %>%
  colnames() %>%
  data.frame(colnames = .) %>%
  filter(str_detect(colnames, ".x|.y")) %>%
  pull(colnames) %>%
  sort()

colnames.x <- colnames.dup %>%
  str_subset(., ".x")

colnames.y <- colnames.dup %>%
  str_subset(., ".y")

new_cols <- colnames.dup %>%
  str_remove(., ".x|.y") %>%
  unique()

## Duplex
phylotable.samples.dup <- phylotable.samples.tmp %>%
  select(all_of(colnames.dup), taxonomy)

## Single
phylotable.samples.sin <- phylotable.samples.tmp %>%
  select(-all_of(colnames.dup), taxonomy)

## Duplex to single
dup.to.sin <- phylotable.samples.dup %>%
  group_by(taxonomy) %>%
  tidyr::pivot_longer(cols = starts_with("MFD"), names_to = "fieldsample_barcode", values_to = "Count") %>%
  ungroup() %>%
  mutate(across(fieldsample_barcode, ~str_remove(., ".x|.y"))) %>%
  group_by(taxonomy, fieldsample_barcode) %>%
  mutate(Count_sum = sum(Count), .keep = "unused") %>%
  ungroup() %>%
  rename(Count = Count_sum) %>%
  distinct() %>%
  group_by(taxonomy) %>%
  tidyr::pivot_wider(names_from = "fieldsample_barcode", values_from = "Count")

## Sort sample IDs
sorted.samples <- mapping.samples %>%
  pull(fieldsample_barcode) %>%
  unique() %>%
  sort()

## Overwrite "OTU" key
otukey <- phylotable.samples.tmp %>%
  select(OTU, Kingdom:Species) %>%
  tidyr::unite(col = "taxonomy", 2:8, sep = ",") %>%
  mutate(across(taxonomy, ~str_replace_all(., "NA", "")))

## Write "OTU" key
data.table::fwrite(otukey, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "N/A",
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_otukey.csv"))


### Make final table
phylotable.samples.final <- phylotable.samples.sin %>%
  left_join(dup.to.sin, by = "taxonomy") %>%
  left_join(otukey, by = 'taxonomy') %>%
  select(OTU, any_of(sorted.samples), taxonomy) %>%
  arrange(., str_rank(OTU, numeric = TRUE)) %>%
  #mutate(across(where(is.numeric), ~tidyr::replace_na(., 0))) %>% # Do this with UNIX sed command
  tidyr::separate_wider_delim(taxonomy, delim = ",",
                              names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                              too_few = "align_start")

# Print to files
data.table::fwrite(phylotable.controls, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "N/A",
                   paste0("release/", format(Sys.time(), "%Y-%m-%d"), "_MFD_arcbac_shallow_release_controls.csv"))

data.table::fwrite(phylotable.samples.final, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "N/A",
                   paste0("release/", format(Sys.time(), "%Y-%m-%d"), "_MFD_arcbac_shallow_release.csv"))

rm(list = ls())
gc()

