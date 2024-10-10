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

## Reproducibility

We include a script [`runme.sh`](./runme.sh) to reproduce the main results of our experiments. This script will run 10 times each necessary execution of the algorithms. Each resulting partition will be stored in the folder `run_artifacts`, the partitions information in `result.csv` and the potential function values in `potential.csv` and `potential2.csv`. In addition it will generate a latex report with all the figures (stored in the `latex_report/img` folder) and tables in the `latex_report` folder. Creating the latex file requires the following dependencies (in addition to the dependencies listed above):

- [latexmk](https://ctan.org/pkg/latexmk/)
- [amsmath](https://ctan.org/pkg/amsmath)
- [amsfonts](https://ctan.org/pkg/amsfonts)
- [multirow](https://ctan.org/pkg/multirow)
- [booktabs](https://ctan.org/pkg/booktabs)
- [section](https://ctan.org/pkg/section)

The `runme.sh` script will reproduce:

- Table 2,3,4,7,8
- Figure 3,4

Table 1 is not reproduced as it is simply fine tuning of parameters. Table 5 and 6 are not reproduced as they do require long computation times. Additionally, Table 5 shows the performances of our algorithms on set system generated from graph which is also the case of Table 3 and 4 that are reproduced.

Figure 1 and 2 are not relevant for reproducibility and Figure 5 is a plot of Table 8 that we reproduce.

We also provide the [`clear_run.sh`](./clear_run.sh) script to clean all data from previous run. It should also be ran once to create the necessary files before running `runme.sh`. Finally we provide the `runme-test.sh` script that will run our algorithm on one instance of the grid set system and generate a report based on that result to check if all the necessary dependencies are installed. If the `runme-test.sh` scripts runs and generates a report with a title and one image. The `runme.sh` script should run without dependencies troubles.

```bash
./clear_run.sh
./runme-test.sh
./clear_run.sh
./runme.sh
```

## Organization

The code is divided in different files which all have their own documentation in the [doc](./doc) folder:

- [main.cpp](./doc/main.md): Script to try out the different variation of the algorithm on different set system
- [partition.cpp](./doc/partition.md): The different variations of the algorithms
- [distances.cpp](./doc/distance.md): Variation of the weight function for the parallel algorithm detailed in the paper
- [classes.cpp](./doc/classes.md): Set systems generation functions
- [visualization](./doc/visualization.md): Scripts to visualize the results obtained

- [main_osm_graphs.cpp](./doc/main_osm_graphs.md): Script to handle graph set systems generated from Open Street Map data
- [main_fb.cpp](./doc/main_fb.md): Script to handle graph set systems generated from graph data such as facebook social circles and ArXiv coauthorship graph
