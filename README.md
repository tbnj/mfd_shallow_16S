# mfd_shallow_16S
The scripts in this repository are part of the [Microflora Danica project](https://github.com/cmc-aau/mfd_wiki/wiki). 
The holds the workflow used to create observational tables of 16S and 18S rRNA gene fragments from raw short-read metagenomes. 

Be advised, that for the metagenomic-derived data, the term "OTU" is only used due to format requirements by ampvis2, and they do not represent classical OTUs. 
The generated profiles can be thought of as taxonomic bins. 

## Scripts
### Bash processing 
`scripts/bcl2fastq.sh` performs demultiplexing using bcl2fastq based on a sample sheet per run. 
The sample sheet file is formatted to output data in the following format:
| Folder | Content |
| --- | --- |
| run folder | A folder of the entire sequencing run. |
| ├── lane folder/ | One folder per flow cell lane |
| │   ├── library folder/ | One folder per library ID |
| │   │   ├── seq files | Two files (forward and reverse) per library ID |


`scripts/trim_rename.sh` performs trimming and quality filtering of the raw data using [fastp](https://doi.org/10.1093/bioinformatics/bty560). Compresses the outputs using [pigz](https://zlib.net/pigz/) as the internal zipping function of fastp is corrupt. 
The output is structured in the following format:
| Folder | Content |
| --- | --- |
| run folder | A folder of the entire sequencing run. |
| ├── lib folder/ | One folder per library ID |
| │   ├── seq files | Two files (forward and reverse) per library ID |


`scripts/classify_shallow_metagenomes.sh` identifies the 16S and 18S using Hidden Markow Models and [nhmmer](https://doi.org/10.1093/bioinformatics/btt403), extracts the reads using awk, grep, sort, uniq and zcat, and performs taxonomic classification using [SINTAX](https://www.drive5.com/usearch/manual/cmd_sintax.html). 
The output is structured in the following format (simpliied to only using results for 16S):
| Folder | Content |
| --- | --- |
| run folder | A folder of the entire sequencing run. |
| ├── lib folder/ | One folder per library |
| │   ├── forward/ |  |
| │   │   ├── arc_LIB-MJXXX-A1_01_forward.hmmout.txt | Output from nhmmer using arc.hmm |
| │   │   ├── bac_LIB-MJXXX-A1_01_forward.hmmout.txt | Output from nhmmer using bac.hmm |
| │   │   ├── arc_bac_LIB-MJXXX-A1_01_forward.fq | Extracted forward reads in fastq format |
| │   │   ├── arc_bac_LIB-MJXXX-A1_01_forward_MFG_ssu_database_NR987_trunc.sintax | Forward SINTAX file |
| │   ├── reverse/ |  |
| │   │   ├── arc_LIB-MJXXX-A1_01_reverse.hmmout.txt | Output from nhmmer using arc.hmm |
| │   │   ├── bac_LIB-MJXXX-A1_01_reverse.hmmout.txt | Output from nhmmer using bac.hmm |
| │   │   ├── arc_bac_LIB-MJXXX-A1_01_reverse.fq | Extracted reverse reads in fastq format |
| │   │   ├── arc_bac_LIB-MJXXX-A1_01_reverse_MFG_ssu_database_NR987_trunc.sintax | Reverse SINTAX file |
| │   ├── tmp/ |  |
| │   │   ├── forward_IDs.txt | Forward read IDs from nhmmer outputs |
| │   │   ├── forward_IDs_unique.txt | Unique forward read IDs |
| │   │   ├── reverse_IDs.txt | Reverse read IDs from nhmmer outputs |
| │   │   ├── reverse_IDs_unique.txt | Unique reverse read IDs |


`scripts/find_hmm.sh` finds all [nhmmer](https://doi.org/10.1093/bioinformatics/btt403) output files and creates symlinks. 


`scripts/find_sintax.sh` finds all [SINTAX](https://www.drive5.com/usearch/manual/cmd_sintax.html) output files and creates symlinks. 

### R processing
`scripts/filter_reads.R` filters the combined pool of 16S rRNA gene reads to within the region between the 8F and 1391R primers for bacteria and 20F and SSU1000ArR for archaea. Splits output based on whether the origin was an MFD sample or an included control. 


`scripts/arcbac_sintax_to_combined_tax.R` creates one big file from all [SINTAX](https://www.drive5.com/usearch/manual/cmd_sintax.html) classifications. Filters the file so that either the forward or the reverse classification file is included if both are present. The read with the deepest classification is used. If the same depth is achieved the forward read is used.  


`scripts/arcbac_combined_tax_to_OTU.R` transforms the outputs of `arcbac_sintax_to_combined_tax.R` to observational tables. The term "OTU" is only used for compatability with [ampvis2](https://kasperskytte.github.io/ampvis2/articles/ampvis2.html). 


To yield the final table N/A is changed to 0 using 
`sed -i 's/N\/A/0/g' mfd_shallow_16S/R/release/DATE__MFD_arcbac_shallow_release.csv`
