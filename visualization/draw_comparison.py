import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd
import numpy as np
import math as m

label_table = ["MP-Mat","minWeight","greedy","sampling","greedy","atOnce","atOnce","atOnce","atOnce","atOnce","atOnce"]
color_table = ["blue", "red","green","orange","green","purple","purple","purple","purple","purple","purple"]
join_table = {'d': 2, 't': 512, 'n':8192}

df = pd.read_csv("experiments_grid.csv", sep=";")
toplot = df.drop(df[(df["t"] != join_table['t']) | (df["d"] != join_table['d']) | (df["ss_type"] != sys.argv[2]) | (df["algo"] == 2)| (df["algo"] == 4)].index).sort_values(sys.argv[1])
if sys.argv[1] == 't':
    toplot = df.drop(df[(df["n"] != join_table['n']) | (df["d"] != join_table['d']) | (df["ss_type"] != sys.argv[2])| (df["algo"] == 2)| (df["algo"] == 4)].index).sort_values(sys.argv[1])
if sys.argv[1] == 'd':
    toplot = df.drop(df[(df["t"] != join_table['t']) | (df["n"] != join_table['n']) | (df["ss_type"] != sys.argv[2])| (df["algo"] == 2)| (df["algo"] == 4)].index).sort_values(sys.argv[1])

fig, ax = plt.subplots()

for key, grp in toplot.groupby('algo'):
    temp = grp.groupby([sys.argv[1]],as_index=False).agg(min=pd.NamedAgg(column="max_crossing", aggfunc="min"), max=pd.NamedAgg(column="max_crossing", aggfunc="max"))
    grp.groupby([sys.argv[1]]).mean("max_crossing").plot(ax=ax,y="max_crossing", label="$\kappa_\mathcal{F}$ " + label_table[int(key)], color=color_table[int(key)])
    ax.fill_between(x=sys.argv[1],y1='min',y2='max',data=temp,color=mpl.colors.to_rgba(color_table[int(key)], 0.15))

    temp = grp.groupby([sys.argv[1]],as_index=False).agg(min=pd.NamedAgg(column="avg_crossing", aggfunc="min"), max=pd.NamedAgg(column="avg_crossing", aggfunc="max"))
    grp.groupby([sys.argv[1]]).mean("avg_crossing").plot(ax=ax,y="avg_crossing", label="$\\bar{\kappa_\mathcal{F}}$ " + label_table[int(key)], color=color_table[int(key)],linestyle='dotted')
    ax.fill_between(x=sys.argv[1],y1='min',y2='max',data=temp,color=mpl.colors.to_rgba(color_table[int(key)], 0.05))


ax.set_xlabel(sys.argv[1], fontsize=25)
ax.set_ylabel("crossing number",fontsize=25)


if sys.argv[1] == 'n':
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[(512**(1-1/2)) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='$512^{1-1/2}$')
    #plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[(2*512**(1-1/2)+ m.log(x,2)) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='2*512^(1-1/2)')
elif sys.argv[1] == 't':
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[(x**(1-1/2) + m.log(8192,2)) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='$t^{1-1/2}$')
    #plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[(x**(1-1/2)+ m.log(8192,2)) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='2*x^(1-1/2)')
elif sys.argv[1] == 'd':
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[(512**(1-1/x)+ m.log(8192,2)) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='$512^{1-1/d}$')
    #plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[(512**(1-1/x)+ m.log(8192,2)) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='2*512^(1-1/d)')

plt.legend(loc='best', fontsize="xx-large")
plt.show()