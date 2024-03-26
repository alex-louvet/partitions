# main.cpp

Script to try out the different variation of the algorithm on different set system

## Arguments

- `-a string` : This is a space separated list of algorithm identifiers to run on the same set system as well as the number of times the given algorithm should be run.
`./.a.out -a '1 1 2 3 5 1'` runs algorithm 1 once, algorithm 2 three times and algorithm 5 once. Algorithm number is between 1 and 11 and are detailes in [partition.cpp](./partition.md)
- `-n int` : Number of points in the set system
- `-t int` : Number of partitions to build
- `-f string` : Type of set system (explained in [classes.cpp](./classes.md)) among `["grid","random_hs","grid_graph","linear_grid","exponential_grid","directional_grid"]`
- `-s` : Save the results. It saves a line in `results.csv` with the parameters of the set system, the number of crossings, average crossing number and the runtime and creates a file with the points and the partition sets
- `-r int` : Seed for random generation 
