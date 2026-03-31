#!/bin/bash
# run_exp1.sh - Experiment 1 automation

GRIDS=("250 100" "500 200" "1000 400")
PARTICLES=(100 10000 1000000 10000000 100000000)

echo "NX,NY,NUM_Points,Approach,Threads,Interp_Time,Mover_Time,Total_Time" > exp1_results.csv

for grid in "${GRIDS[@]}"; do
    set -- $grid
    NX=$1
    NY=$2
    
    for NP in "${PARTICLES[@]}"; do
        for APPROACH in 1 2; do
            # Use sed to modify main.cpp temporarily
            sed -i "s/NX       = .*/NX       = $NX;/" main.cpp
            sed -i "s/NY       = .*/NY       = $NY;/" main.cpp
            sed -i "s/NUM_Points = .*/NUM_Points = $NP;/" main.cpp
            sed -i "s/approach = .*/approach = $APPROACH;/" main.cpp
            sed -i "s/num_threads = .*/num_threads = 1;/" main.cpp
            
            g++ main.cpp init.cpp utils.cpp -lm -fopenmp -O2 -o main
            OUTPUT=$(./main | grep "CSV:")
            echo "$OUTPUT" | sed 's/CSV: //' >> exp1_results.csv
            echo "Done: NX=$NX, NY=$NY, NP=$NP, Approach=$APPROACH"
        done
    done
done

echo "Experiment 1 complete! Results saved to exp1_results.csv"