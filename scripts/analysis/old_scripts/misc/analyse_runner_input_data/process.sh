#!/bin/bash
input=input.data
rm -f process.out
rm -f energy.out atom_info.out
grep "\benergy\b" $input >> energy.out
grep "\batom\b" $input >> atom_info.out
python3.7 analysis.py $input >> process.out
grep "begin" $input | wc -l
grep "atom" $input | wc -l
eog 01_energy_histogram.png
