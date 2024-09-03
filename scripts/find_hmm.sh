#!/usr/bin/env bash
#### Script to search for and the creation of symbolic links to nhmmer output ####

# paths
DIR_IN=('/mfd_shallow_16S/data')
DIR_OUT_F=('/mfd_shallow_16S/data/hmm_forward_out/')
DIR_OUT_R=('/mfd_shallow_16S/data/hmm_reverse_out/')

# find samples - add min and max depth  
find $DIR_IN -type f -name 'arc_*_forward.hmmout.txt' -exec ln -s '{}' $DIR_OUT_F ';'
find $DIR_IN -type f -name 'bac_*_forward.hmmout.txt' -exec ln -s '{}' $DIR_OUT_F ';'
find $DIR_IN -type f -name 'euk_*_forward.hmmout.txt' -exec ln -s '{}' $DIR_OUT_F ';'
find $DIR_IN -type f -name 'arc_*_reverse.hmmout.txt' -exec ln -s '{}' $DIR_OUT_R ';'
find $DIR_IN -type f -name 'bac_*_reverse.hmmout.txt' -exec ln -s '{}' $DIR_OUT_F ';'
find $DIR_IN -type f -name 'euk_*_reverse.hmmout.txt' -exec ln -s '{}' $DIR_OUT_F ';'
