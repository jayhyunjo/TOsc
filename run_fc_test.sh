#!/bin/bash
cd /home/jjo/Work/MicroBooNE/ToAll_basic_nospeedup
for ROW in 35 47 52 56; do
  rm -rf test_fc_row$ROW
  for j in $(seq 0 18); do
    ./scan_3p1_exact -j $j -n 19 -t 7000 -m 500 -k 50 -gdm2 60 -gamp 60 \
        -row $ROW -crn 1 -fc 1 -o test_fc_row$ROW > /dev/null 2>&1 &
  done
  wait
  echo "FC row $ROW done"
done
date "+FC TEST DONE %H:%M %Z"
