import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np

label_table = ["","min","rate","sampling","deque","atOnce","atOnce","atOnce","atOnce","atOnce","atOnce"]
color_table = ["", "red","green","orange","blue","purple","purple","purple","purple","purple","purple"]
join_table = {'d': 2, 't': 512, 'n':8192}

df = pd.read_csv("temp.csv", sep=";")
toplot = df.drop(df[(df["t"] != join_table['t']) | (df["d"] != join_table['d']) | (df["ss_type"] != int(sys.argv[2])) ].index).sort_values(sys.argv[1])
if sys.argv[1] == 't':
    toplot = df.drop(df[(df["n"] != join_table['n']) | (df["d"] != join_table['d']) | (df["ss_type"] != int(sys.argv[2]))].index).sort_values(sys.argv[1])
if sys.argv[1] == 'd':
    toplot = df.drop(df[(df["t"] != join_table['t']) | (df["n"] != join_table['n']) | (df["ss_type"] != int(sys.argv[2]))].index).sort_values(sys.argv[1])

fig, ax = plt.subplots()

for key, grp in toplot.groupby(['algo']):
    ax = grp.groupby([sys.argv[1]]).mean().plot(ax=ax, kind='line', y="time", label=label_table[int(key)], color=color_table[int(key)])

ax.set_xlabel(sys.argv[1], fontsize=25)
ax.set_ylabel("runtime", fontsize=25)

plt.legend(loc='best', fontsize="xx-large")
plt.show()
