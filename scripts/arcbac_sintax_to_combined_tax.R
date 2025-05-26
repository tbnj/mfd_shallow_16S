### Setup env
#library(tidyverse)
library(dplyr)
library(stringr)

setwd("/mfd_shallow_16S")

### Import data

##Import best hit from hmms
arc_reads <- readr::read_delim("output/DATE_ARC_reads.txt", col_names = c("Read_name"), delim = ",")
bac_reads <- readr::read_delim("output/DATE_BAC_reads.txt", col_names = c("Read_name"), delim = ",")

arcbac_reads <- rbind(arc_reads, bac_reads) %>%
  pull(Read_name)

## Import SINTAX data
forward_files <- list.files("sintax_forward_out", pattern = "arc_bac", full.names = T)
reverse_files <- list.files("sintax_reverse_out", pattern = "arc_bac", full.names = T)

## Read all forward files
forward <- do.call(rbind, Map(cbind, lapply(forward_files, data.table::fread, sep = "\t", header = FALSE, fill = TRUE), 
                              V5 = basename(forward_files))) %>%
  rename(Read_name = 1, Tax_string = 2, Strand = 3, Taxonomy = 4, File_name = 5) %>%
  mutate(Read_name = str_remove(Read_name, "\\s.*")) %>%
  select(-Strand) %>%
  filter(Read_name %in% arcbac_reads) %>%
  mutate(across(Taxonomy, ~na_if(., ""))) %>%
  filter(!is.na(Taxonomy))

## Read all reverse files
reverse <- do.call(rbind, Map(cbind, lapply(reverse_files, data.table::fread, sep = "\t", header = FALSE, fill = TRUE), 
                              V5 = basename(reverse_files))) %>%
  rename(Read_name = 1, Tax_string = 2, Strand = 3, Taxonomy = 4, File_name = 5) %>%
  mutate(Read_name = str_remove(Read_name, "\\s.*")) %>%
  select(-Strand) %>%
  filter(Read_name %in% arcbac_reads) %>%
  mutate(across(Taxonomy, ~na_if(., ""))) %>%
  filter(!is.na(Taxonomy))

## Write temporary output
data.table::fwrite(forward, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_forward_arcbac_reads_sintax.csv"))

data.table::fwrite(reverse, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_reverse_arcbac_reads_sintax.csv"))


## Identify paied reads found both in forward and reverse
pairs <- Reduce(intersect, list(forward$Read_name, reverse$Read_name))


### Create subsets 

## Build list
lst_df <- lst(forward, reverse)

## Select taxonomic information
lst_df.tax <- lst_df %>%
  lapply(., select, c(1,3))

## Filter for paired reads
tmp_pairs <- lst_df.tax %>%
  bind_rows(.id = 'Source') %>%
  filter(Read_name %in% pairs)

## Filter for non-paired reads
tmp_unique <- lst_df.tax %>%
  bind_rows(.id = 'Source') %>%
  filter(!Read_name %in% pairs)

## Transform paired to long format
tmp_pairs.long <- tmp_pairs %>%
  group_by(Read_name) %>%
  distinct(Taxonomy, .keep_all = TRUE) %>%
  mutate(duplicate = n()) %>%
  ungroup()

## Unique reads forward
unique_forward <- tmp_unique %>%
  filter(Source == "forward")

## Unique reads reverse
unique_reverse <- tmp_unique %>%
  filter(Source == "reverse")

## Identify pairs with taxonomic agreement
same_name <- tmp_pairs.long %>%
  filter(!duplicate == 2) %>% 
  select(-duplicate)

## Identify pairs with taxonomic disagreement - select best (if same depth slice chooses first entry = forward)
diff_name <- tmp_pairs.long %>%
  filter(!duplicate == 1) %>% 
  select(-duplicate) %>%
  #mutate(rowsum = rowSums(is.na(.)))
  mutate(rowsum = str_count(Taxonomy, ":")) %>%
  group_by(Read_name) %>%
  slice(which.max(rowsum)) %>%
  select(-rowsum) %>%
  ungroup() %>%
  mutate(Source = "combined")

# Test
# tmp <- tmp_pairs.long %>%
#   filter(!duplicate == 1) %>% 
#   select(-duplicate) %>%
#   mutate(rowsum = str_count(Taxonomy, ":")) %>%
#   group_by(Read_name, rowsum) %>%
#   mutate(duplicate = n()) %>%
#   ungroup() %>%
#   filter(duplicate == 2) %>%
#   group_by(Read_name) %>%
#   slice(which.max(rowsum))


### Combine outputs and minimise
df_combined <- bind_rows(same_name, diff_name, unique_forward, unique_reverse) %>%
  tidyr::separate_wider_delim(Taxonomy, delim = ",", 
                              names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), too_few = "align_start") %>%
  mutate(across(Kingdom:Species, ~str_remove(., ".:"))) %>%
  select(-Source)

### Write to disk
out_path <- paste0(getwd(), "/sintax_combined_out/")

data.table::fwrite(df_combined, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0(out_path, format(Sys.time(), "%Y-%m-%d"), "_arcbac_sintax_combined.csv"))

rm(list = ls())
gc()




