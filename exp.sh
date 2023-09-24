#!/bin/bash
DLIST=(2 10 100)
for ((p = 1024; p < 17000; p*=16)); do
    for ((t = 16; t < 65; t*=2)); do
        for d in "${DLIST[@]}"; do
            ./a.out 1 $p $t 1 $d
        done
    done
done
