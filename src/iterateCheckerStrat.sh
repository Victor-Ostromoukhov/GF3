#!/usr/bin/env zsh


for ind1 in {1..10}
do
  for ind2 in {1..10}
  do
     sem -j+0 ./checkerStrat2D -i $ind1 -j $ind2 --seed1min 1 --seed1max 10 --seed2min 1 --seed2max 10
  done
done
sem --wait
