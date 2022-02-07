#!/bin/bash

# intention: submit script for MD simulations of different thermalized graphene slabs for PES development

cd heatup_350K
  nohup ../heatup_300K/md_tian2.x md_tian.inp &
  wait

cd ../heatup_400K
  nohup ../heatup_300K/md_tian2.x md_tian.inp &
  wait

cd ../heatup_450K
  nohup ../heatup_300K/md_tian2.x md_tian.inp &
  wait

cd ../heatup_500K
  nohup ../heatup_300K/md_tian2.x md_tian.inp &
  wait

cd ../heatup_550K
nohup ../heatup_300K/md_tian2.x md_tian.inp &
