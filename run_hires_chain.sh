#!/bin/bash
# Waits for the running 60x60 corrected scan to finish, then launches the
# high-resolution overnight run:
#   amplitude 120 bins (8%/bin, was 16.6%) -> quantization below the MC floor
#   dm2       100 bins (25/decade, was 15)  -> resolves the 1-10 eV2 structure
#   1000 SM toys (median noise 1.0%, was 1.4%)
#   7000 p-value toys per hypothesis (publication count), common random numbers
cd /home/jjo/Work/MicroBooNE/ToAll_basic_nospeedup
while pgrep -f "scan_3p1_exact [-]j" > /dev/null; do sleep 60; done
sleep 30
rm -rf output_3p1_hires logs_3p1_hires
mkdir -p logs_3p1_hires
date "+HIRES START %H:%M %Z"
for j in $(seq 0 18); do
  ./scan_3p1_exact -j $j -n 19 -t 7000 -m 1000 -k 50 -gdm2 100 -gamp 120 -crn 1 \
      -o output_3p1_hires > logs_3p1_hires/job_$j.log 2>&1 &
done
wait
date "+HIRES DONE %H:%M %Z"
echo "files: $(ls output_3p1_hires | wc -l)"
