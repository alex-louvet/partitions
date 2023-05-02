#generates pair object of a given ground set X
def generate_pairs(X):
    import classes
    res = []
    for i in range(len(X)-1):
        for j in range(i+1,len(X)):
            res.append(classes.Edge.frompointlist((X[i],X[j])))
    return res

#Returns whether and edge separates a set (i.e. one element in, one out)
def intersects(e,s):
    first_in_s = False
    second_in_s = False
    for x in s.points:
        if compare(x,e.points[0]):
            first_in_s = True
        if compare(x,e.points[1]):
            second_in_s = True
    return first_in_s != second_in_s

#check equality of vectors
def compare(x,y):
    if len(x) != len(y):
        return False
    for i in range(len(x)):
        if x[i] != y[i]:
            return False
    return True

def sample_edge(omega,radius):
    import random as r

    #determine highest non-empty weight bucket
    start = 0
    while len(omega[start]) == 0:
        start += 1

    #compute summed weight of buckets
    index = start
    pow = 1
    sum = 0
    while index < len(omega) and index - start < radius:
        sum += len(omega[index])*pow
        index += 1
        pow /= 2

    #randomly pick from weight
    index = start
    random = r.random()*sum
    partial_sum = 0
    pow = 1
    while partial_sum < random:
        partial_sum += len(omega[index])*pow
        index += 1
        pow /= 2

    return r.choice(omega[index-1])

def sample_set(pi,radius):
    import random as r
    import math as m

    index = m.floor(max(0,len(pi)-radius))

    #compute summed weight of bucket
    pow = 1
    sum = 0
    while index < len(pi):
        sum += len(pi[index])*pow
        index += 1
        pow *= 2

    #randomly pick from weight
    index = min(0,len(pi))
    random = r.random()*sum
    partial_sum = 0
    pow = 1
    while partial_sum < random:
        partial_sum += len(pi[index])*pow
        index += 1
        pow *= 2

    return r.choice(pi[index-1])

def build_matching(X, S, a,b,gamma, radius_e, radius_s, warmup,threads):
    import classes
    M=[]
    i = 1
    while len(X) >= 4:
        print('iteration',i)
        Mp = match_half(X,S,a*len(X)**gamma+b, radius_e, radius_s, warmup,threads)
        M = M+Mp
        newX = []
        for x in X:
            test = False
            for e in Mp:
                if compare(x,e.points[0]) or compare(x,e.points[1]):
                    test = True
            if not test:
                newX.append(x)
        X = newX
        print('')
        i += 1
    if len(X) > 1:
        M.append(classes.Edge.frompointlist(X))
    return M

def match_half(X,S,kappa, radius_e, radius_s, warmup, threads):
    import math
    import random as r
    from tqdm import tqdm
    from joblib import Parallel, delayed

    n = len(X)
    m = len(S)
    E = n*(n-1)/2

    p = min(106*n/kappa**2*math.log(E*n/4),1)
    q = min(39*n/kappa**2*math.log(m*n/4),1)

    omega = [generate_pairs(X),[]]
    pi = [S,[]]

    res = []

    for i in tqdm(range(warmup + math.ceil(n/4))):
        #sample one edge and set according to weights
        ei = sample_edge(omega, radius_e)
        si = sample_set(pi, radius_s)

        #sample sets and edges uniformly w/ proba p and q
        Ei = []
        Si = []
        for j in range(len(omega)):
            newL = []
            for e in omega[j]:
                if r.random() <= p:
                    Ei.append(e)
                else:
                    newL.append(e)
            omega[j] = newL

        for j in range(len(pi)):
            newL = []
            for e in pi[j]:
                if r.random() <= q:
                    Si.append(e)
                else:
                    newL.append(e)
            pi[j] = newL


        #update weight depending on intersection with ei,si
        temp = Parallel(n_jobs=threads)(delayed(intersects)(e,si) for e in Ei)
        for j,e in enumerate(Ei):
            if temp[j]:
                e.increase()
            while len(omega) - 1 < e.intersection:
                omega.append([])
            omega[e.intersection].append(e)
        temp = Parallel(n_jobs=threads)(delayed(intersects)(ei,s) for s in Si)
        for j,s in enumerate(Si):
            if temp[j]:
                s.increase()
            while len(pi) - 1 < s.intersection:
                pi.append([])
            pi[s.intersection].append(s)

        #remove ei and adjacent edges from weights
        if i >= warmup:
            res.append(ei)
            for j in range(len(omega)):
                newL = []
                for e in omega[j]:
                    if not compare(e.points[0], ei.points[0]) and not compare(e.points[1], ei.points[1]) and not compare(e.points[0], ei.points[1]) and not compare(e.points[1], ei.points[0]):
                        newL.append(e)
                omega[j] = newL

    return res

def check_intersect(matching, S, save, name, warmup, threads):
    from joblib import Parallel, delayed
    import csv
    import shutil

    intersect_table = [0]*len(matching)
    for i,e in enumerate(matching):
        intersect_para = Parallel(n_jobs=threads)(delayed(intersects)(e,s) for s in S)
        intersect_table[i] = sum(intersect_para)

    if save:
        with open(name + '_' + str(warmup) + '.csv', 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(intersect_table)
            shutil.move(name + '_' + str(warmup) + '.csv', 'results/')
    return intersect_table


def just_train_weights(X,S,kappa, radius_e, radius_s, iterations, threads, start_weight):
    import math
    import random as r
    from tqdm import tqdm
    from joblib import Parallel, delayed

    n = len(X)
    m = len(S)
    E = n*(n-1)/2

    p = min(106*n/kappa**2*math.log(E*n/4),1)
    q = min(39*n/kappa**2*math.log(m*n/4),1)

    omega = [generate_pairs(X),[]]
    pi = start_weight

    for _ in tqdm(range(iterations)):
        #sample one edge and set according to weights
        ei = sample_edge(omega, radius_e)
        si = sample_set(pi, radius_s)

        #sample sets and edges uniformly w/ proba p and q
        Ei = []
        Si = []
        for j in range(len(omega)):
            newL = []
            for e in omega[j]:
                if r.random() <= p:
                    Ei.append(e)
                else:
                    newL.append(e)
            omega[j] = newL

        for j in range(len(pi)):
            newL = []
            for e in pi[j]:
                if r.random() <= q:
                    Si.append(e)
                else:
                    newL.append(e)
            pi[j] = newL

        #update weight depending on intersection with ei,si
        temp = Parallel(n_jobs=threads)(delayed(intersects)(e,si) for e in Ei)
        for j,e in enumerate(Ei):
            if temp[j]:
                e.increase()
            while len(omega) - 1 < e.intersection:
                omega.append([])
            omega[e.intersection].append(e)
        temp = Parallel(n_jobs=threads)(delayed(intersects)(ei,s) for s in Si)
        for j,s in enumerate(Si):
            if temp[j]:
                s.increase()
            while len(pi) - 1 < s.intersection:
                pi.append([])
            pi[s.intersection].append(s)

    return omega,pi
