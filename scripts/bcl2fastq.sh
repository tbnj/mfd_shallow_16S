#!/usr/bin/env bash
#### Run bcl2fastq for individual sequencing runs ####

## CAPITALISATION = CHANGE
RUNNAME=RUN_JOURNAL
RUNFOLDER=/PATH_TO_RAW_DIR/${RUNNAME}/
OUTDIR=/PATH_TO_OUT_DIR/${RUNNAME}
SAMPLESHEET=/PATH_TO_SAMPLE_SHEET/DATE_JOURNAL.csv
THREADS_PROC=20
THREADS_LOAD=echo "$((THREADS/5))"
THREADS_WRITE=echo "$((THREADS/5))"

bcl2fastq \
--runfolder-dir $RUNFOLDER \
--output-dir $OUTDIR \
--ignore-missing-bcls \
--ignore-missing-filter \
--ignore-missing-positions \
--mask-short-adapter-reads 0 \
--sample-sheet $SAMPLESHEET \
--loading-threads $THREADS_LOAD \
--processing-threads $THREADS_PROC \
--writing-threads $THREADS_WRITE \
&> /tmp/${RUNNAME}.log && \
mv /tmp/${RUNNAME}.log /$OUTDIR/${RUNNAME}.log
