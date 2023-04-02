#!/bin/bash
LOG_DIR='logs'
mkdir $LOG_DIR
LOGFILE=$LOG_DIR/`date +"%Y%m%d-%H%M"`.log

# INPUT VARS
INPUT=individual_primer/*

for PROGRAM in '_boxplots.py' '_PCoA.py'
do
    echo "[LOG] `date`: Started $PROGRAM !" 2>&1 | tee -a $LOGFILE
    (time (\
        python src/$PROGRAM \
            $INPUT \
            --exclude-primer 'off-target')
    ) 2>&1 | tee -a $LOGFILE
    echo "[LOG] `date`: Finished $PROGRAM !" 2>&1 | tee -a $LOGFILE
done