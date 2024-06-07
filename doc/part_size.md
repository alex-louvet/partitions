# part_size.cpp

This file regroups the different functions that build an integer table representing the number of elements in each part of the partition.

- 10 : `vector<int> equi_distro(int n, int t)` builds a table of size `t` with n/t in each cell. It represents a partition where all parts have the same size.

- `vector<int> reduce_size(int n, int t, int min)` builds a table with where half of the points will be in partition of size n/t, then 1/4 in partition of size n/(2t) etc... until parts of size `min` that will contain all the remaining points.

- 11 : `vector<int> interval_size(int n, int t, int min, int max)` builds a table with where half of the points will be in partition of size `max`, then 1/4 in partition of size `max`/2 etc... until parts of size `min` that will contain all the remaining points.

- 12 : `vector<int> same_number_different_size_3(int n, int t)` builds a table with one third of the parts have size 3n/2t, one third have size n/t and the remaining third have size n/2t.

- 13 : `vector<int> same_number_different_size_linear(int n, int t, float mineps)` builds a table with `t` parts of size evenly spaced from (1+`mineps`)n/t to (1-`mineps`)n/t