# partition.cpp

Contains the different algorithms presented in the paper and additional ones

- 1 : min (in the paper). Picks the point of minimum weight at each iteration and updates the weight of sets at each iteration.
- 2 : rate (in the paper). Picks any point with weight smaller than thy detailed in the paper at each iteration and updates the weight of sets at each iteration.
- 3 : sampling (not included in the paper because it gave bad results) -> In this variation, we sample sets at each round and we compute point weights only on these sets. We then compare it to the (scaled) rate described in the paper. 
- 4 : Dequeue algorithm (in the paper). Picks any point with weight smaller than thy detailed in the paper at each iteration and updates the weight of sets only if no point with small weight can be found.
- 5-10 : Parallel algorithm (in the paper). Build partition at once based on a weights computed after each partition is built. We tried different [distance.cpp](./distance.md) functions to give weight to points to build a given partition.
