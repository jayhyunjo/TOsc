#!/bin/bash
# Final production run: confirmed Nature inputs, production lattice, capping
# toys, 10k pseudos, 1000 SM toys, and the global-ensemble seed fix.
cd /home/jjo/Work/MicroBooNE/ToAll_basic_nospeedup
rm -rf output_3p1_final2 logs_3p1_final2
mkdir -p logs_3p1_final2
date "+FINAL2 START %H:%M %Z" > final2.log
for j in $(seq 0 18); do
  ./scan_3p1_exact -j $j -n 19 -t 10000 -m 1000 -k 50 -gdm2 80 -gamp 60 \
      -toyclip 1 -o output_3p1_final2 > logs_3p1_final2/job_$j.log 2>&1 &
done
wait
date "+FINAL2 DONE %H:%M %Z" >> final2.log
/home/jjo/root-build/bin/root -l -b -q 'plot_3p1_exact.C("output_3p1_final2", 80, 60)' > extract_final2.txt 2>&1
mv contour_3p1_exact.root contour_3p1_final2.root 2>/dev/null
mv contour_3p1_exact.png  contour_3p1_final2.png  2>/dev/null
date "+FINAL2 EXTRACTION DONE %H:%M %Z" >> final2.log
