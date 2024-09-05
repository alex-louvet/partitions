# A Greedy Algorithm for Low-Crossing Partitions for General Set Systems

Monika Csikos, Alexandre Louvet, Nabil H. Mustafa

This is the code for the experiments related to the paper "A Greedy Algorithm for Low-Crossing Partitions for General Set Systems" accepted to ALENEX25. We provide C++ implementations of the code detailed in the paper as well as some python visualization scripts.

## Dependencies

### C++ dependencies

- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [OpenMP](https://www.openmp.org/)

### Python dependencies

- [Matplotlib](https://matplotlib.org/)
- [Pandas](https://pandas.pydata.org/)
- [Numpy](https://numpy.org/)
- [Scipy](https://scipy.org/)

## Organization

The code is divided in different files which all have their own documentation in the [doc](./doc) folder:

- [main.cpp](./doc/main.md): Script to try out the different variation of the algorithm on different set system
- [partition.cpp](./doc/partition.md): The different variations of the algorithms
- [distances.cpp](./doc/distance.md): Variation of the weight function for the parallel algorithm detailed in the paper
- [classes.cpp](./doc/classes.md): Set systems generation functions
- [visualization](./doc/visualization.md): Scripts to visualize the results obtained

- [main_osm_graphs.cpp](./doc/main_osm_graphs.md): Script to handle graph set systems generated from Open Street Map data
- [main_fb.cpp](./doc/main_fb.md): Script to handle graph set systems generated from graph data such as facebook social circles and ArXiv coauthorship graph
