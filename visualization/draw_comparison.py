import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np

label_table = ["","min","","rate","","","","","dequeue","","","","","","","","","","parallel"]
color_table = ["", "red","", "green","","","","","blue","","","","","","","","","","purple"]
join_table = {'d': 2, 't': 512, 'n':8192}

df = pd.read_csv("results.csv", sep=";")

toplot = df.drop(df[(df["t"] != join_table['t']) | (df["d"] != join_table['d']) | (df["ss_type"] != sys.argv[2]) ].index).sort_values(sys.argv[1])
if sys.argv[1] == 't':
    toplot = df.drop(df[(df["n"] != join_table['n']) | (df["d"] != join_table['d']) | (df["ss_type"] != sys.argv[2])].index).sort_values(sys.argv[1])
if sys.argv[1] == 'd':
    toplot = df.drop(df[(df["t"] != join_table['t']) | (df["n"] != join_table['n']) | (df["ss_type"] != sys.argv[2])].index).sort_values(sys.argv[1])

print(toplot)

fig, ax = plt.subplots()

for key, grp in toplot.groupby(['algo']):
    print(key)
    ax = grp.groupby([sys.argv[1]]).mean().plot(ax=ax, kind='line', y="max_crossing", label=label_table[int(key)], color=color_table[int(key)])
    ax = grp.groupby([sys.argv[1]]).mean().plot(ax=ax, kind='line', y="avg_crossing", label='avg ' + label_table[int(key)], color=color_table[int(key)],linestyle="dotted")

if sys.argv[1] == 'n':
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[512**(1-1/2) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='512^(1-1/2)')
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[2*512**(1-1/2) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='2*512^(1-1/2)')
elif sys.argv[1] == 't':
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[2*x**(1-1/2) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='x^(1-1/2)')
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[x**(1-1/2) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='2*x^(1-1/2)')
elif sys.argv[1] == 'd':
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[2*512**(1-1/x) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='512^(1-1/d)')
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[512**(1-1/x) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='2*512^(1-1/d)')
plt.legend(loc='best')
plt.show()
