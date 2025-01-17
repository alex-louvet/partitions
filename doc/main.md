# main.cpp

Script to try out the different variation of the algorithm on different set system

Compile with

```bash
g++ -fopenmp main.cpp
```

## Arguments

- `-a string` : This is a space separated list of algorithm identifiers to run on the same set system as well as the number of times the given algorithm should be run.
`./.a.out -a '1 1 2 3 5 1'` runs algorithm 1 once, algorithm 2 three times and algorithm 5 once. Algorithm number is between 1 and 10 and are detailed in [partition.cpp](./partition.md)
- `-n int` : Number of points in the set system. Default: `1024`
- `-t int` : Number of partitions to build. Default: `16`
- `-f string` : Type of set system (explained in [classes.cpp](./classes.md)) among `["grid","grid_centers","random_hs","grid_graph","linear_grid","exponential_grid","directed_grid","random","projective_plane","ERGGraph","power_law","concentric_circles,"file"]`. Default: `grid`
- `-d int`: Dimension in which the points of $X$ live or distance to build the adjecency set system from the graph in case of graph set systems. Explained in more details in [classes.cpp](./classes.md). Default `2`
- `-c int`: Number of centers to generate when giving option `-f grid_centers
- `-s` : Save the results. It saves a line in `results.csv` with the parameters of the set system, the number of crossings, average crossing number and the runtime and creates a file with the points and the partition sets
- `-r int` : Seed for random generation. Default: `time(NULL)`
- `-k int` : Number of rounds to run the different sample algorithms as described in the paper. Default `log(m*n)`
- `-i string` : Name of the file containing the set system when giving option `-f file`. Specification for the files can be found in [classes.cpp](./classes.md)
- `-e` : Also compute the maximum value of $\lvert |F|/|X| - |A\cap F|/|A|\rvert$ overall $F$ of $\mathcal{F}$ with $A$ a random sample of $X$ and $A$ composed of one element of each part of the computed partition. This allows to compare the approximation performances of approximation obtained via random sample and  thy obtained from a partition

for instance if compiled with the command above

```bash
./a.out -a '1 1 3 4' -n 4096 -t 45 -d 7 -f grid -s
```

will run algorithm 1 (min) once and algorithm 3 (sampling) four times on a set system consisting of 4096 random points and 7-dimensional axis-orthogonal halfspaces evenly spaced. It will create 44 partitions of 91 points and one partition of 92 points. It will save the crossing number obtained in `results.csv` and will create 5 files containing the points of the set system and the partition sets as incidence vectors obtained for each iteration of one of the algorithms.