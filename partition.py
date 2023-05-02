def first_algo(X,S,t,radius_e,radius_s,warmup,threads):
    import matching_algo as match
    import math
    import random as r
    kept = X
    D = len(X[0])
    m = len(S)
    partitiontable = {}

    # create a hashtable with key the elts of X and value a list with only themselves
    for x in X:
        partitiontable[tuple(x)] = [x]

    # Iterated matchings
    while len(kept) > t:
        print('partition:',len(kept))
        temp = match.build_matching(kept,S,10**(1/D),math.log(m),1-(1/D), radius_e, radius_s, warmup,threads)
        kept=[]
        for e in temp:
            if r.random() < 0.5:
                kept.append(e.points[0])
                partitiontable[tuple(e.points[0])] = partitiontable[tuple(e.points[1])] + partitiontable[tuple(e.points[0])]
                del partitiontable[tuple(e.points[1])]
            else:
                kept.append(e.points[1])
                partitiontable[tuple(e.points[1])] = partitiontable[tuple(e.points[1])] + partitiontable[tuple(e.points[0])]
                del partitiontable[tuple(e.points[0])]

    return [partitiontable[i] for i in partitiontable]

def intersects(p,s):
    import matching_algo as match
    inside = False
    for y in s.points:
        if match.compare(p[0],y):
            inside = True
            break
    for x in p[1:]:
        x_in = False
        for y in s.points:
            if match.compare(x,y):
                if not inside:
                    return True
                else:
                    x_in = True
                    break
        if not x_in and inside:
            return True
    return False

def check_intersect(partitiontable, S, save, name, warmup, threads,variant):
    from joblib import Parallel, delayed
    import csv
    import shutil

    intersect_table = [0]*len(S)
    for i,s in enumerate(S):
        intersect_para = Parallel(n_jobs=threads)(delayed(intersects)(p,s) for p in partitiontable)
        intersect_table[i] = sum(intersect_para)

    if save:
        with open(name + '_' + str(warmup) + '_partition' + variant +'.csv', 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(intersect_table)
            shutil.move(name + '_' + str(warmup) + '_partition' + variant + '.csv', 'results/')
    return intersect_table

def check_intersect_against_grid(partitiontable,save,name,warmup,threads,variant):
    import csv
    import shutil
    import math as m

    n = sum([len(part) for part in partitiontable])

    intersect_table = [0]*(2*m.ceil(m.sqrt(n)))
    for part in partitiontable:
        maxx = max([x[0] for x in part])
        minx = min([x[0] for x in part])
        maxy = max([x[1] for x in part])
        miny = min([x[1] for x in part])
        for i in range(m.ceil(minx*m.sqrt(n)),m.ceil(maxx*m.sqrt(n))):
            intersect_table[i] += 1
        for i in range(m.ceil(miny*m.sqrt(n)),m.ceil(maxy*m.sqrt(n))):
            intersect_table[i+m.ceil(m.sqrt(n))] += 1

    if save:
        with open(name + '_' + str(warmup) + '_partition' + variant +'.csv', 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(intersect_table)
            shutil.move(name + '_' + str(warmup) + '_partition' + variant + '.csv', 'results/')
    return intersect_table


def second_algo(X,S,t,radius_e,radius_s,iterations,threads):
    import math
    import matching_algo as match
    from tqdm import tqdm
    res = []
    D = len(X[0])
    m = len(S)
    n = len(X)
    weight = [S,[]]
    for k in range(t):
        print('partition:',k+1)
        partition = []
        omega,pi = match.just_train_weights(X,S,10**(1/D)*len(X)**(1-1/D)+math.log(m), radius_e, radius_s, iterations, threads, weight)
        for _ in tqdm(range(n//t-1)):
            crossing = []
            not_crossing = S
            index = 0
            
            #pick an edge with highest weight and connected
            while len(omega[index]) == 0 :
                index += 1
            if len(partition) == 0:
                e = omega[index].pop(0)
                partition.append(e.points[0])
                partition.append(e.points[1])
            else:
                found = False
                while not found:
                    if index > len(omega):
                        raise ValueError('Went through omega without finding connected edge')
                    for j,edge in enumerate(omega[index]):
                        if found:
                            break
                        for x in partition:
                            if match.compare(x,edge.points[0]):
                                found = True
                                e = omega[index].pop(j)
                                partition.append(edge.points[1])
                                break
                            elif match.compare(x,edge.points[1]):
                                found = True
                                e = omega[index].pop(j)
                                partition.append(edge.points[0])
                                break
                    index += 1

                #Delete edges linking newly added point and points already in the partition
                for L in omega:
                    temp = []
                    for i,e in enumerate(L):
                        if match.compare(partition[-1],e.points[0]):
                            for x in partition:
                                if match.compare(x,e.points[1]):
                                    temp.append(i)

                        if match.compare(partition[-1],e.points[1]):
                            for x in partition:
                                if match.compare(x,e.points[0]):
                                    temp.append(i)
                    for x in temp[::-1]:
                        del L[x]

            #Select sets that picked edge crosses
            temp = []
            for s in not_crossing:
                if match.intersects(e,s):
                    crossing.append(s)
                else:
                    temp.append(s)
            not_crossing = temp

            #increase weight of edges crossed by lines crossing picked edge
            for i,L in enumerate(omega[1:]):
                newL = []
                for e in L:
                    crosses = False
                    for s in crossing:
                        if match.intersects(e,s):
                            e.decrease()
                            omega[i].append(e)
                            crosses = True
                            break
                    if not crosses:
                        newL.append(e)
                omega[i+1] = newL

        #Build initial weights of next iteration
        weight = []
        for s in crossing:
            s.algo2inc()
            while len(weight) <= s.algo2numhascrossed:
                weight.append([])
            weight[s.algo2numhascrossed].append(s)
        for s in not_crossing:
            while len(weight) <= s.algo2numhascrossed:
                weight.append([])
            weight[s.algo2numhascrossed].append(s)

        temp = []
        for x in X:
            keep = True
            for y in partition:
                if match.compare(x,y):
                    keep = False
                    break
            if keep:
                temp.append(x)
        X = temp
        res.append(partition)
    return res

def intersects_return_weight(e,s):
    import matching_algo as match
    first_in_s = False
    second_in_s = False
    for x in s.points:
        if match.compare(x,e.points[0]):
            first_in_s = True
        if match.compare(x,e.points[1]):
            second_in_s = True
    return (first_in_s != second_in_s)*s.algo3weight

def update_weights_algo3(e,S):
    res = 0
    for s in S:
        res += intersects_return_weight(e,s)
    return res

def compute_intersection_matrix_row(e,S):
    import matching_algo as match
    import numpy as np
    res = np.zeros(len(S),dtype=bool)
    for s in S:
        res[s.id] = match.intersects(e,s)
    return res

def third_algo(X,S,t,threads):
    from joblib import Parallel, delayed
    import matching_algo as match
    from tqdm import tqdm
    import time
    import numpy as np
    print('Algo3')
    res = []
    n = len(X)
    edges_alive = match.generate_pairs(X)
    M = np.matrix(Parallel(n_jobs=threads)(delayed(compute_intersection_matrix_row)(edges_alive[i],S) for i in tqdm(range(len(edges_alive)))),dtype=bool)
    for k in range(t):
        print('')
        print('partition',k+1)
        partition = []
        setsWeight = np.array([s.algo3weight for s in S])
        for i,e in enumerate(edges_alive):
            e.algo3weight = np.sum(np.dot(M[e.id], setsWeight))
        crossing = []
        not_crossing = S
        for i in tqdm(range(n//t - 1)):
            #select edge of minimum weight
            if partition == []:
                select_from = [(e,2) for e in edges_alive]
            else:
                select_from = []
                for e in edges_alive:
                    for p in partition:
                        if match.compare(e.points[0],p):
                            select_from.append((e,1))
                        elif match.compare(e.points[1],p):
                            select_from.append((e,0))
            min = select_from[0]
            for i,e in enumerate(select_from):
                if e[0].algo3weight < min[0].algo3weight:
                    min = e
            if partition == []:
                partition.append(min[0].points[0])
                partition.append(min[0].points[1])
            else:
                partition.append(min[0].points[min[1]])

            # remove edges connecting 2 points in the partition
            temp = []
            for i,e in enumerate(edges_alive):
                if match.compare(partition[-1],e.points[0]):
                    for x in partition:
                        if match.compare(x,e.points[1]):
                            temp.append(i)

                if match.compare(partition[-1],e.points[1]):
                    for x in partition:
                        if match.compare(x,e.points[0]):
                            temp.append(i)
            for x in temp[::-1]:
                del edges_alive[x]

            # Update weight of edges with crossing hs
            still_not_crossing = []
            for s in not_crossing:
                if M[min[0].id,s.id]:
                    crossing.append(s)
                    for e in edges_alive:
                        if M[e.id,s.id]:
                            e.algo3weight -= s.algo3weight
                    s.algo3inc()
                else:
                    still_not_crossing.append(s)
            not_crossing = still_not_crossing

        #remove edges connected to points in the partition
        temp = []
        for i,e in enumerate(edges_alive):
            for x in partition:
                if match.compare(x,e.points[1]) or match.compare(x,e.points[0]):
                    temp.append(i)
        for x in temp[::-1]:
            del edges_alive[x]

        res.append(partition)

    print('')
    return res

def draw_partition(partition, name):
    import scipy
    import matplotlib.pyplot as plt
    import math as m
    n = sum([len(part) for part in partition])
    for k in range(m.ceil(m.sqrt(n))+1):
        plt.plot([k/m.ceil(m.sqrt(n)),k/m.ceil(m.sqrt(n))],[0,1], color='lightgrey')
    for k in range(m.ceil(m.sqrt(n))+1):
        plt.plot([0,1],[k/m.ceil(m.sqrt(n)),k/m.ceil(m.sqrt(n))], color='lightgrey')
    
    if len(partition[0][0]) != 2:
        raise ValueError('Wrong dimension to plot: required 2, got ', len(partition[0][0]))
    for i,p in enumerate(partition):
        ch = scipy.spatial.ConvexHull(p)
        plt.plot([x[0] for x in p], [x[1] for x in p], 'o', label='partition' + str(i+1))
        for simplex in ch.simplices:
            plt.plot([p[simplex[0]][0],p[simplex[1]][0]], [p[simplex[0]][1], p[simplex[1]][1]], 'k-')
    # plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.15), ncols=4, fancybox=True, shadow=True)
    plt.savefig(name)

def list_grid_intersects(e,n):
    import math as m
    res = []
    if e.points[0][0] < e.points[1][0]:
        minx = e.points[0][0]
        maxx = e.points[1][0]
    else:
        minx = e.points[1][0]
        maxx = e.points[0][0]
    if e.points[0][1] < e.points[1][1]:
        miny = e.points[0][1]
        maxy = e.points[1][1]
    else:
        miny = e.points[1][1]
        maxy = e.points[0][1]
    for i in range(m.ceil(minx*m.sqrt(n)),m.ceil(maxx*m.sqrt(n))):
        res.append(i)
    for i in range(m.ceil(miny*m.sqrt(n)),m.ceil(maxy*m.sqrt(n))):
        res.append(i+m.ceil(m.sqrt(n)))
    return res

def intersects_grid(e,k,n):
    import math as m
    if k < m.ceil(m.sqrt(n)):
        if e.points[0][0] < e.points[1][0]:
            minx = e.points[0][0]
            maxx = e.points[1][0]
        else:
            minx = e.points[1][0]
            maxx = e.points[0][0]
        return minx*m.ceil(m.sqrt(n)) <= k and k <= maxx*m.ceil(m.sqrt(n))
    else:
        if e.points[0][1] < e.points[1][1]:
            miny = e.points[0][1]
            maxy = e.points[1][1]
        else:
            miny = e.points[1][1]
            maxy = e.points[0][1]
        return miny*m.ceil(m.sqrt(n)) <= k-m.ceil(m.sqrt(n)) and k <= maxy*m.ceil(m.sqrt(n))


def third_algo_grid(X,S,t,threads):
    from joblib import Parallel, delayed
    import matching_algo as match
    from tqdm import tqdm
    import numpy as np
    print('Algo3')
    res = []
    n = len(X)
    edges_alive = match.generate_pairs(X)
    for k in range(t):
        print('')
        print('partition',k+1)
        partition = []
        for i,e in enumerate(edges_alive):
            e.algo3weight = len(list_grid_intersects(e,n))
        crossing = []
        not_crossing = S
        for i in tqdm(range(n//t - 1)):
            #select edge of minimum weight
            if partition == []:
                select_from = [(e,2) for e in edges_alive]
            else:
                select_from = []
                for e in edges_alive:
                    for p in partition:
                        if match.compare(e.points[0],p):
                            select_from.append((e,1))
                        elif match.compare(e.points[1],p):
                            select_from.append((e,0))
            min = select_from[0]
            for i,e in enumerate(select_from):
                if e[0].algo3weight < min[0].algo3weight:
                    min = e
            if partition == []:
                partition.append(min[0].points[0])
                partition.append(min[0].points[1])
            else:
                partition.append(min[0].points[min[1]])

            # remove edges connecting 2 points in the partition
            temp = []
            for i,e in enumerate(edges_alive):
                if match.compare(partition[-1],e.points[0]):
                    for x in partition:
                        if match.compare(x,e.points[1]):
                            temp.append(i)

                if match.compare(partition[-1],e.points[1]):
                    for x in partition:
                        if match.compare(x,e.points[0]):
                            temp.append(i)
            for x in temp[::-1]:
                del edges_alive[x]

            # Update weight of edges with crossing hs
            still_not_crossing = []
            for s in not_crossing:
                if M[min[0].id,s.id]:
                    crossing.append(s)
                    for e in edges_alive:
                        if M[e.id,s.id]:
                            e.algo3weight -= s.algo3weight
                    s.algo3inc()
                else:
                    still_not_crossing.append(s)
            not_crossing = still_not_crossing

        #remove edges connected to points in the partition
        temp = []
        for i,e in enumerate(edges_alive):
            for x in partition:
                if match.compare(x,e.points[1]) or match.compare(x,e.points[0]):
                    temp.append(i)
        for x in temp[::-1]:
            del edges_alive[x]

        res.append(partition)

    print('')
    return res
