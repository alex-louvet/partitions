import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str)
parser.add_argument("-t", type=int)

args=parser.parse_args()
file= args.file
t = args.t

algo1 = [0]*(t+1)

with open(file, newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in spamreader:
        for x in row:
            algo1[int(x)] += 1

plt.bar([i for i in range(t+1)],algo1,color='red', label='algo 1')
plt.show()
