# Simplicial partition construction point by point

This is the code for the experiments related to the paper [Simplicial partition construction point by point](link to come). We provide C++ implementations of the code detailed in the paper as well as some python visualization scripts.

## Dependencies

We detail the dependencies for the code to be run

### C++ dependencies

- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [OpenMP](https://www.openmp.org/)

### Python dependencies

- [Matplotlib](https://matplotlib.org/)
- [Pandas](https://pandas.pydata.org/)
- [Numpy](https://numpy.org/)
- [Scipy](https://scipy.org/)

## Organisation

The code is divided in different files which all have their own documentation in the [doc](./doc) folder:

- [main.cpp](./doc/main.md): Script to try out the different variation of the algorithm on different set system
- [partition.cpp](./doc/partition.md): The different variations of the algorithms
- [distances.cpp](./doc/distance.md): Variation of the weight function for the parallel algorithm detailed in the paper
- [classes.cpp](./doc/classes.md): Set systems generation functions
- [visualization](./doc/visualization.md): Scripts to visualize the results obtained
