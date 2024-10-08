# distances.cpp

This files regroups the different distance function explored for the parallel algorithm of the paper. They correspond to arguments 5 to 11 int argument `-a` of [main](./main.md)

- 5 : L1 distance
- 6 : L2 distance
- 7 : We sample sets and add 1 to the weights of points which edge with the root of the partition intersects this set
- 8 : We maintain weights on points and sets and each round sample one point and one set according to their respective distribution. We then increase the weight of point from which the sample set intersect their edge with the root and decrease the weight of the sets which intersect the edge between the root and the sampled point. This is inspired from work on [Low crossing number matchings](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2021.28)
- 9 : Same as 7 but instead of adding 1 to the weight we add the weight of the set obtained from previous partition building
- 10-13 : Variation described in the paper with different [part_size](./part_size.md)
