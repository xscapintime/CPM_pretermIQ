#!/bin/bash

seeds=($(seq 1 20 1000))

for s in ${seeds[@]}
do
    echo $s
    #python 2a.corr_iteration.py $s
    python test_iter.py $s

done
