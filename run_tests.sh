#!/bin/bash

# Define arrays for each parameter
THREAD_COUNTS=(4 8 16)
FRAME_NUMBERS=(1 2 4 8 16 32 64 128)
INPUT_FILES=("in_5K.fluid" "in_15K.fluid" "in_35K.fluid" "in_100K.fluid" "in_300K.fluid" "in_500K.fluid")

# Use nested loops to run the app with all combinations of parameters
for FILE in "${INPUT_FILES[@]}"; do
    for THREAD in "${THREAD_COUNTS[@]}"; do
        for FRAME in "${FRAME_NUMBERS[@]}"; do
            for RUN in {1..5}; do
                echo "Running with THREAD=$THREAD, FRAME=$FRAME, FILE=$FILE"
                ./fluidanimate $THREAD $FRAME $FILE
            done 
        done
    done
done

echo "All combinations completed!"

