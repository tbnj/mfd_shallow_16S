#!/usr/bin/env bash
#### Script to search for and the creation of symbolic links to SINTAX output ####

# Paths
DIR_IN=('/mfd_shallow_16S/data')
DIR_OUT_F=('/mfd_shallow_16S/data/sintax_forward_out/')
DIR_OUT_R=('/mfd_shallow_16S/data/sintax_reverse_out/')

# Find results for ARC and BAC
find $DIR_IN -type f -name 'arc_bac*_forward_MFD_ssu_database_NR987_trunc.sintax' -exec ln -s '{}' $DIR_OUT_F ';'
find $DIR_IN -type f -name 'arc_bac*_reverse_MFD_ssu_database_NR987_trunc.sintax' -exec ln -s '{}' $DIR_OUT_R ';'

# Find results for EUK
find $DIR_IN -type f -name 'euk*_forward_MFD_ssu_database_NR987_trunc.sintax' -exec ln -s '{}' $DIR_OUT_F ';'
find $DIR_IN -type f -name 'euk*_reverse_MFD_ssu_database_NR987_trunc.sintax' -exec ln -s '{}' $DIR_OUT_R ';'