import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd
import numpy as np

label_table = ["","min","rate","sampling","deque","atOnce","atOnce","atOnce","atOnce","atOnce","atOnce"]
color_table = ["", "red","green","orange","blue","purple","purple","purple","purple","purple","purple"]
join_table = {'d': 2, 't': 512, 'n':8192}

df = pd.read_csv("results_17_04.csv", sep=";")
toplot = df.drop(df[(df["t"] != join_table['t']) | (df["d"] != join_table['d']) | (df["ss_type"] != sys.argv[2]) | (df["algo"] > 4)].index).sort_values(sys.argv[1])
if sys.argv[1] == 't':
    toplot = df.drop(df[(df["n"] != join_table['n']) | (df["d"] != join_table['d']) | (df["ss_type"] != sys.argv[2]) | (df["algo"] > 4)].index).sort_values(sys.argv[1])
if sys.argv[1] == 'd':
    toplot = df.drop(df[(df["t"] != join_table['t']) | (df["n"] != join_table['n']) | (df["ss_type"] != sys.argv[2]) | (df["algo"] > 4)].index).sort_values(sys.argv[1])

fig, ax = plt.subplots()

for key, grp in toplot.groupby(['algo']):
    temp = grp.groupby([sys.argv[1]],as_index=False).agg(min=pd.NamedAgg(column="rate_violation", aggfunc="min"), max=pd.NamedAgg(column="rate_violation", aggfunc="max"))
    grp.groupby([sys.argv[1]]).mean("rate_violation").plot(ax=ax,y="rate_violation", label=label_table[int(key)], color=color_table[int(key)])
    ax.fill_between(x=sys.argv[1],y1='min',y2='max',data=temp,color=mpl.colors.to_rgba(color_table[int(key)], 0.15))


ax.set_xlabel(sys.argv[1], fontsize=25)
ax.set_ylabel("rate_violation",fontsize=25)

plt.legend(loc='best', fontsize="xx-large")
plt.show()