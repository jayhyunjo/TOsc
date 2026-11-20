#!/bin/bash
# Runs while the user is away. Waits for the row matrix, extends both candidate
# configurations to two more dm2 rows, then writes a single summary table.
cd /home/jjo/Work/MicroBooNE/ToAll_basic_nospeedup
LOG=away.log
note(){ date "+%H:%M:%S $1" >> $LOG; }

note "waiting for row matrix"
while [ "$(ps -C scan_3p1_exact -o pid= | wc -l)" -gt 0 ]; do sleep 120; done
sleep 20
note "row matrix finished"

run(){ TAG=$1; ROW=$2; shift 2
  rm -rf test_${TAG}_r${ROW}
  for j in 0 1 2; do
    ./scan_3p1_exact -j $j -n 3 -t 10000 -m 1000 -k 50 -gdm2 80 -gamp 60 \
        -toyclip 1 -row $ROW "$@" -o test_${TAG}_r${ROW} > /dev/null 2>&1 &
  done
  wait
  note "$TAG row $ROW done"
}

# extend both candidates to dm2 ~ 0.30 (row 30) and ~ 29.9 (row 72)
for ROW in 30 72; do run fcboth $ROW -fc 1; done
for ROW in 30 72; do run asimov $ROW -prof 1; done

note "writing summary"
/home/jjo/root-build/bin/root -l -b -q compare_rows.C > away_summary.txt 2>&1
note "AWAY QUEUE COMPLETE"
