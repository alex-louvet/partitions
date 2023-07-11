#!/bin/bash
./a.out 1 $1 &
./a.out 2 $1 &
./a.out 3 $1 &
./a.out 4 $1 &
./a.out 5 $1 &
./a.out 6 $1 &
wait
