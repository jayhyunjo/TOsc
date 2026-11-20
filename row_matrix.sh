#!/bin/bash
# Row tests at three dm2 values against Nature: control (40), kink (47), high (61).
# Two remaining hypothesis choices, both with BOTH branches (the data-matching config):
#   A: per-toy profiling (fc)            B: Asimov-profiled fixed hypothesis
cd /home/jjo/Work/MicroBooNE/ToAll_basic_nospeedup
run(){ TAG=$1; shift; ROW=$1; shift
  rm -rf test_${TAG}_r${ROW}
  for j in 0 1 2; do
    ./scan_3p1_exact -j $j -n 3 -t 10000 -m 1000 -k 50 -gdm2 80 -gamp 60 \
        -toyclip 1 -row $ROW "$@" -o test_${TAG}_r${ROW} > /dev/null 2>&1 &
  done
  wait
  echo "$TAG row $ROW done" >> rowmatrix.log
}
for ROW in 47 61 40; do run fcboth $ROW -fc 1; done
for ROW in 47 61 40; do run asimov $ROW -prof 1; done
echo "ROWMATRIX DONE" >> rowmatrix.log
