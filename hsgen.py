#generate n random points of dimension d
def generate_points(n,d):
    import random as r
    import numpy as np
    res = []
    for _ in range(n):
        temp = np.zeros(d)
        for i in range(d):
            temp[i] = r.random()
        res.append(tuple(temp))
    return res

def generate_points_on_circle(n):
    import random as r
    import math as m
    res = []
    for _ in range(n):
        theta = r.random()*2*m.pi
        res.append((0.5+(0.5+r.random()*-1*0.05)*m.cos(theta),0.5+(0.5+r.random()*-1*0.05)*m.sin(theta)))
    return res


#compute a halfspace of passing through the points of L
def compute_halfspace(points):
    import numpy as np
    centroid = np.mean(points, axis=0)
    _, _, vh = np.linalg.svd(points - centroid, full_matrices=False)
    hyperplane = vh[-1]
    orientation = np.zeros(len(points[0]))
    orientation[-1] = 1
    if np.dot(hyperplane,orientation) < 0:
        return centroid, -1*vh[-1]
    return centroid, vh[-1]

def check_side(centroid, equation, point):
    import numpy as np
    return np.dot(equation, point - centroid) > 0

def halfspace(points, centroid, equation):
    import classes
    res1 = []
    res2 = []
    for x in points:
        if check_side(centroid,equation,x):
            res1.append(x)
        else:
            res2.append(x)
    return classes.Set.frompointlist(res1), classes.Set.frompointlist(res2)

def generate_random_halfspaces(points,number):
    d = len(points[0])
    res = []
    for _ in range(number//2):
        temp = generate_random_hs(points,d)
        res.append(temp[0])
        res.append(temp[1])
    return res

def generate_random_halfspaces_parallel(points,number, threads):
    from joblib import Parallel, delayed
    d = len(points[0])
    temp = Parallel(n_jobs=threads)(delayed(generate_random_hs)(points,d) for _ in range(number//2))
    res = []
    for x in temp:
        res.append(x[0])
        res.append(x[1])
    return res

def generate_random_hs(points,d):
    import random as r
    temp = r.sample(points,r.randint(1,d))
    c,h = compute_halfspace(temp)
    hs1,hs2 = halfspace(points,c,h)
    return hs1,hs2

def generate_hs_grid(points):
    import math as m
    import classes
    sqr = m.ceil(m.sqrt(len(points)))
    d = len(points[0])
    temp = [[] for _ in range(d*sqr)]
    for x in points:
        for k in range(d):
            for i in range(m.ceil(x[k]*sqr),sqr):
                temp[k*sqr + i].append(x)
    return [classes.Set.frompointlist(s) for s in temp]

# DANGER ZONE
def list_all_halfspaces(points):
    res = []
    import itertools
    d = len(points[0])
    temp = list(itertools.combinations(points,d))
    for comb in temp:
        c,h = compute_halfspace(comb)
        temp = halfspace(points,c,h)
        res.append(temp)
    return res
