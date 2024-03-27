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

### Method

`void buildAdjacency(bool pts)` calls `buildAdjacency` for each set and creates `belongs_to` and `not_in` for each point if `pts = true`

### Instances

`SetSystem Grid(int n, int d)` builds a grid sets system on `n` points and $dn^{1/d}$ sets  evenly spaced and facing towards $(1\ldots 1)$

`SetSystem RandomHyperplanes(int n, int d, int m)` build a random halfspaces set system that are each defined by the hyperplane passing through a random subsets of `d` points and oriented towards $(1\ldots 1)$

`SetSystem GridGraph(int n, int d)` builds `n*n` points on the integer grid and s one set is added to the set system for each point containing points with $L_1$ distance `d` to the considered point

`SetSystem LinearGrid(int n, int d)` similar to the grid set system but the ith set orthogonal to a canonical vector is duplicated i times

`SetSystem ExponentialGrid(int n, int d)` similar to the grid set system but the ith set orthogonal to a canonical vector is duplicated $e^i$ times

`SetSystem DirectionalGrid(int n, int d)` similar to the grid set system the sets orthogonal to $1 0 \ldots 0$ are duplicated $n^{1/d}$ times
