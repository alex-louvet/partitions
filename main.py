import matching_algo as match
import hsgen
import math
import time as t
import pickle
import shutil
import argparse
import partition
import classes

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

# Matching demo
# matching = match.build_matching(X,S,10**(1/D),math.log(m),1-(1/D), radius_e, radius_s, warmup,threads)
#
# intersect_table = match.check_intersect(matching, S, save, name, warmup, threads)
# print(intersect_table)

# Partition demo
# var1=t.time()
# part1 = partition.first_algo(X,S,4,radius_e,radius_s,warmup,threads)
# var1 = t.time() - var1
# var2 = t.time()
# part2 = partition.second_algo(X,S,4,radius_e,radius_s,warmup,threads)
# var2 = t.time()-var2
# var3 = t.time()
part3 = partition.third_algo(X,S,32,threads)
# var3 = t.time() - var3
# var4 = t.time()
# part4 = partition.fourth_algo(X,S,16,1/16,threads)
# var4 = t.time() - var4

# if save:
#     with open(name +'_partitions.pkl','wb') as f:
#         pickle.dump((part1,part2,part3),f)
#     shutil.move(name + '_partitions.pkl','results/')



# intersect_table1 = partition.check_intersect(part1, S, save, name, warmup, threads, '1')
# intersect_table2 = partition.check_intersect(part2, S, save, name, warmup, threads, '2')
intersect_table3 = partition.check_intersect_against_grid(part3, save, name, warmup, threads, '3_grid')
intersect_table3 = partition.check_intersect(part3,S, save, name, warmup, threads, '3_train')
# intersect_table4 = partition.check_intersect(part4, S, save, name, warmup, threads, '4')

# partition.draw_partition(part1)
# partition.draw_partition(part2)
partition.draw_partition(part3, name+'.png')
# partition.draw_partition(part4)


# D = {}
# var = t.time()
# E = match.generate_pairs(X)
# print(t.time()-var)
# var = t.time()
# temp = Parallel(n_jobs=threads)(delayed(partition.compute_intersection_matrix_row)(e,S) for e in E)
# print(t.time()-var)
# var = t.time()
# for i,e in enumerate(E):
#     D[e.points] = temp[i]
# print(t.time()-var)
# var = t.time()
# new = D[E[0].points][S[0].points]
# print(t.time()-var)
# var = t.time()
# new = match.intersects(E[0],S[0])
# print(t.time()-var)
# print([s.id for s in S])
