#!/bin/bash

#EXP ON d
./a.out 1 8192 32 1 2 &
./a.out 1 8192 32 1 3 &
./a.out 1 8192 32 1 4 &
./a.out 1 8192 32 1 5 &
./a.out 1 8192 32 1 10 &
./a.out 2 8192 32 1 2 &
./a.out 2 8192 32 1 3 &
./a.out 2 8192 32 1 4 &
./a.out 2 8192 32 1 5 &
./a.out 2 8192 32 1 10 &
./a.out 3 8192 32 1 2 &
./a.out 3 8192 32 1 3 &
./a.out 3 8192 32 1 4 &
./a.out 3 8192 32 1 5 &
./a.out 3 8192 32 1 10 &

# EXP on n
./a.out 1 2048 32 1 2 &
./a.out 1 4096 32 1 2 &
./a.out 1 8192 32 1 2 &
./a.out 1 16384 32 1 2 &
./a.out 1 32768 32 1 2 &
./a.out 2 2048 32 1 2 &
./a.out 2 4096 32 1 2 &
./a.out 2 8192 32 1 2 &
./a.out 2 16384 32 1 2 &
./a.out 2 32768 32 1 2 &
./a.out 3 2048 32 1 2 &
./a.out 3 4096 32 1 2 &
./a.out 3 8192 32 1 2 &
./a.out 3 16384 32 1 2 &
./a.out 3 32768 32 1 2 &

# EXP on t

./a.out 1 8192 16 2 2 &
./a.out 1 8192 32 2 2 &
./a.out 1 8192 64 2 2 &
./a.out 1 8192 128 2 2 &
./a.out 2 8192 16 2 2 &
./a.out 2 8192 32 2 2 &
./a.out 2 8192 64 2 2 &
./a.out 2 8192 128 2 2 &
./a.out 3 8192 16 2 2 &
./a.out 3 8192 32 2 2 &
./a.out 3 8192 64 2 2 &
./a.out 3 8192 128 2 2 &
