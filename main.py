import matching_algo as match
import hsgen
import math
import time as t
import pickle
import shutil
import argparse
import partition
import classes
from joblib import delayed,Parallel

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str)
parser.add_argument("-n", type=int)
parser.add_argument("-d", type=int)
parser.add_argument("-m", type=int)
parser.add_argument("-t", "--threads", type=int)
parser.add_argument("-w", "--warmup", type=int)
parser.add_argument("-s", "--save", action='store_true')
parser.add_argument("-p", "--partitions", type=int)

args=parser.parse_args()
file = args.file
N = args.n
D = args.d
m = args.m
threads = args.threads
warmup = args.warmup
save = args.save
partitions = args.partitions

if file:
    with open(file,'rb') as f:
        X,S = pickle.load(f)
    tempS = []
    for s in S:
        tempS.append(classes.Set.frompointlist(s.points))
    S = tempS
    N = len(X)
    D = len(X[0])
    m = len(S)
    if not threads:
        threads = 1
    if not warmup:
        warmup = 0
    if not partitions:
        partitions = math.floor(math.sqrt(N))
    name = file[:-4]
else:
    if not N:
        N = 500
    if not D:
        D = 4
    if not m:
        m = 10000
    if not threads:
        threads = 1
    if not warmup:
        warmup = 0
    if not partitions:
        partitions = math.floor(math.sqrt(N))

    # Set system generation
    var = t.time()
    print('Generating points:')
    # X = hsgen.generate_points_on_circle(N)
    X = hsgen.generate_points(N,D)
    print('time:', t.time()-var)
    print('')
    print('Generating halfspaces:')
    #S = hsgen.generate_random_halfspaces(X,m)
    S = hsgen.generate_hs_grid(X)
    m = len(S)
    print('time:', t.time()-var)
    print('')
    name = str(math.floor(t.time())) + '_' + str(N) + '_' + str(D)+ '_' + str(m)+ '_' + str(partitions)
    if save:
        with open(name + '.pkl','wb') as f:
            pickle.dump((X,S),f)
        shutil.move(name + '.pkl','ss_examples/')

radius_e = math.log(N*(N+1)/2)
radius_s = math.log(m*N/4)

def run_algo(X,S,M,version):
    match version:
        case '2':
            return partition.third_algo_random(X,S,partitions,M)
        case '3':
            return partition.third_algo_random_weights(X,S,partitions,M)
        case '4':
            return partition.third_algo_bfs(X,S,partitions,M)
        case '5':
            return partition.third_algo_bfs_min(X,S,partitions,M)
        case '6':
            return partition.third_algo_min_bfs(X,S,partitions,M)
        case '1':
            return partition.third_algo(X,S,partitions,M)

M = partition.parallel_intersection_matrix(X,S,threads) 
min,uniform,weighted,bfs,bfs_min,min_bfs = Parallel(n_jobs=threads)(delayed(run_algo)(X,S,M,version) for version in ['1', '2', '3','4','5','6'])

def histogram(table):
    res = []
    for x in table:
        if len(res) < x + 1:
            res += [0] * (x + 1 - len(res))
        res[x] += 1
    return res

stats1 = histogram(partition.check_intersect_against_grid(min[0], save, name, warmup, threads, 'min'))
stats2 = histogram(partition.check_intersect_against_grid(uniform[0], save, name, warmup, threads, 'uniform'))
stats3 = histogram(partition.check_intersect_against_grid(weighted[0], save, name, warmup, threads, 'weighted'))
stats4 = histogram(partition.check_intersect_against_grid(bfs[0], save, name, warmup, threads, 'bfs'))
stats5 = histogram(partition.check_intersect_against_grid(bfs_min[0], save, name, warmup, threads, 'bfs_min'))
stats6 = histogram(partition.check_intersect_against_grid(min_bfs[0], save, name, warmup, threads, 'min_bfs'))

partition.draw_partition(min[0], name+'_min.png',stats1,min[1] , 'Min algo')
partition.draw_partition(uniform[0], name+'_uniform.png',stats2, uniform[1] ,'Uniform w/ rate')
partition.draw_partition(weighted[0], name+'_weighted.png',stats3, weighted[1] ,'Weights w/ rate')
partition.draw_partition(bfs[0], name+'_bfs.png',stats4, bfs[1] ,'Arbitrary BFS')
partition.draw_partition(bfs_min[0], name+'_bfsmin.png',stats5, bfs_min[1] ,'BFS + min')
partition.draw_partition(min_bfs[0], name+'_minbfs.png',stats6, min_bfs[1] ,'Min + BFS')
