# classes.cpp

This is the file referring to set system generation, it contains all the object defined for the purpose of the experiments.

## `class Point`

### Attributes

`vector<float> coordinates` (size `d`) contains the coordinates of the point

`vector<float> belongs_to` (size `m`) list of sets, represented by their number, the point belongs to

`vector<float> not_in` (size `m`) list of sets, represented by their number, that the point doesn't belongs to

### Constructor

`Point(int dimension)` construction with random coordinates

`Point(vector<float> coords)` construction with given coordinates

## `class Edge`

### Attribute

`int points[2]` tuple of points froming the Edge

`int weight` weight of the edge

### Constructor

`Edge(int index1, int index2)` construction by giving the two index of points in the set system

## `class Set`

### Attributes

`vector<bool> points` cine of the incidence matrix corresponding to the set

`int weight` weight of the set

`vector<int> points_indices` list of points, represented by their number, that belong to the set

`vector<int> complement_indices` list of points, represented by their number, that don't belong to the set

### Constructor

`Set(vector<bool> indices)` built by giving the adjacency list

`Set(int size)` builds an empty set of a set system with `n` points

### Method

`void buildAdjacency()` builds `points_indices` and `complement_indices` from the incidence vector `points`

## `class SetSystem`

### Attributes

`vector<Point> points` list of points

`vector<Set> sets` list of sets

### Constructor

`SetSystem(int n,int d,int m)` build set system with `n` random points in dimension `d` and `m` empty sets

`SetSystem(vector<Point> p,vector<Set> s)` builds set system with a given list of points and list of sets

`SetSystem()` builds empty set system

`SetSystem(filename string)` builds the set system from the file `filename`. The file must contain first the points of the set system with a point per line where coordinates of the point are flllowed by a comma `,` (even the last one). A line containing the word `sets` separates the points from the adjacency vector of the sets. After that each line represents a set with a comma after each element of the adjacency vector. Each set must be represented by a vector of size $n$. The file [setSystem.example](./setSystem.example) illustrates how to format the file.

### Method

`void buildAdjacency(bool pts)` calls `buildAdjacency` for each set and creates `belongs_to` and `not_in` for each point if `pts = true`

### Instances

`SetSystem Grid(int n, int d)` builds a grid sets system on `n` points and $dn^{1/d}$ sets  evenly spaced and facing towards $(1\ldots 1)$

`SetSystem GridWithCenters(int n, int d, int c)` builds a grid sets system on `n` points and $dn^{1/d}$ sets  evenly spaced and facing towards $(1\ldots 1)$ but points are generated following a gaussian distribution around `c` different centers generated uniformly. Each center gets the same number of points generated from the distribution centered on them

`SetSystem RandomHyperplanes(int n, int d, int m)` build a random halfspaces set system that are each defined by the hyperplane passing through a random subsets of `d` points and oriented towards $(1\ldots 1)$

`SetSystem GridGraph(int n, int d)` builds `n*n` points on the integer grid and s one set is added to the set system for each point containing points with $L_1$ distance `d` to the considered point

`SetSystem LinearGrid(int n, int d)` similar to the grid set system but the ith set orthogonal to a canonical vector is duplicated i times

`SetSystem ExponentialGrid(int n, int d)` similar to the grid set system but the ith set orthogonal to a canonical vector is duplicated $e^i$ times

`SetSystem DirectionalGrid(int n, int d)` similar to the grid set system the sets orthogonal to $1 0 \ldots 0$ are duplicated $n^{1/d}$ times

`SetSystem Random(int n, int d, int m, float p)` generates a random set system with $n$ points in dimension $d$ and $m$ sets where each point has probbility $p$ to belong to each set.

`SetSystem ProjectivePlane(int n)` generates a projective plane of order n where n is a prime. That is a ground set of size $n^2+n+1$ and sets with the properties of [projective planes](https://en.wikipedia.org/wiki/Projective_plane). The code to generate them is adapted from [Salmelu's work](https://github.com/Salmelu/ProjectivePlane).

`SetSystem ERGGraph(int n, int d, float p)` generates a set system from an [Erdos–Rényi–Gilbert graph](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model) $G$ with probability p for an edge to exist and where for each element $x$ of $X$, there exists a set $F = \{y \in X \vert D_G(x,y) \le d\}$ where $D_G$ is the distance in $G$.

`SetSystem PowerLaw(int n, int d, float p)` generates a [power law graph](https://mathweb.ucsd.edu/~fan/power.pdf) $G$ with the degree of each node random with probability of degree $k$ proportional to $1/k^{\beta}$ and a random matching of nodes nodes to mtach their randomly decided degree. The code to generate the graph is adapted from [Coudert, Csikos, Ducoffe and Viennot's work](https://gitlab.inria.fr/viennot/graph-vcdim).

`SetSystem ConcentricCircles(int n, int m)` generates 4 concentric circles of radius 0.6, 0.7, 0.8 and 0.9 centered on the origin and containing each 1/4th of the `n` elements of $X$ . Then the sets are determined by disks centered at a random elements of $X$ and with a radius determined by drawing another element of $X$. 