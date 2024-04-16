# main_osm_graphs.cpp

Script to try out the different variation of the algorithm on different set system

Compile with

```bash
g++ -fopenmp main_osm_graphs.cpp
```

## Arguments

- `-a string` : This is a space separated list of algorithm identifiers to run on the same set system as well as the number of times the given algorithm should be run.
`./.a.out -a '1 1 2 3 5 1'` runs algorithm 1 once, algorithm 2 three times and algorithm 5 once. Algorithm number is between 1 and 10 and are detailed in [partition.cpp](./partition.md)
- `-n string` : Name of the set system. Default `dakota`
- `-t int` : Number of partitions to build. Default: `16`
- `-d int` : Distance to build the adjacency set system from the graph. Default `1`
- `-s` : Save the results. It saves a line in `results.csv` with the parameters of the set system, the number of crossings, average crossing number and the runtime and creates a file with the points and the partition sets
- `-r int` : Seed for random generation. Default: `time(NULL)`
- `-k int` : Number of rounds to run the different sample algorithms as described in the paper. Default `log(m*n)`
- `-e string` : Name of the file containing the edges for the graph
- `-g string` : Name of the file containing the gps coordinates of the graph

for instance if compiled with the command above

```bash
./a.out -a '1 1 3 4' -n paris -t 16 -d 1 -g data/paris-geo.txt -e data/paris-edges.txt -s
```

will run algorithm 1 (min) once and algorithm 3 (sampling) four times on a set system consisting of the graph represented in `data/paris-edges.txt` and `data/paris-geo.txt` with distances `1`, i.e. sets are simply edges. It will save the crossing number obtained in `results.csv` and will create 5 files containing the points of the set system and the partition sets as incidence vectors obtained for each iteration of one of the algorithms.