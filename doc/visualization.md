# Visualization scripts

These scripts can be found in the [visualization](./visualization) folder.

`draw_comparison.py [n,t,d] string set_system` draw the graph for the variation of the crossing number depending on either `n`, `t` or `d` for the given type of `set_system`. By default, it gets the data from `results.csv`

![result image of draw_comparison.py](../img/grid.png)

`draw_runtime.py [n,t,d] string set_system` draw the graph for the variation of the runtime depending on either `n`, `t` or `d` for the given type of `set_system`. By default, it gets the data from `results.csv`

![result image of draw_runtime.py](../img/gridrt.png)

`draw_partition.py string file (y)` draws the partition from a partition `file`. Add the argument `y` if you don't want the lightgrey grid to be drawn

![result image of draw_partition.py](../img/gridr.png)

`draw_partition_seq.py string file` draws the partitions from a partition `file` one after the other

`draw_partition_osm.py string file string edges (y)` draws the partition for a graph set system obtained from osm data from a partition `file`. It will draw the edges stores in the `edges` file. Add the argyment `y` if you want the outline of each partition to be drawn as well

![result image of draw_partition_osm.py](../img/osm.png)

`approx_graphs.py file` draws the graph of error factor from file `file`

![result image of approx_graphs.py](../img/approx_t.png)

`draw_activity.py file1 file2` draw bar graphs from a list of number of violations representain the number of violation averaged over the runs with the same parameters

![result image of draw_activity.py](../img/cob_16s.png)

`draw_violation.py file` draw number of violation in a graph (replaced in the paper by bar plot generated with `draw_activity.py`)