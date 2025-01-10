import math as m
import sys
from logging import raiseExceptions

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

name = sys.argv[1]

dlist = [2, 3, 4, 5, 10]
tlist = [16, 32, 64, 128, 256]

var = sys.argv[2]
if var not in ["d", "t"]:
    raise ValueError("Second arg must be t or d")

label_table = [
    "",
    "MinWeight",
    "GreedyPotential",
    "sampling",
    "GreedyPotential",
    "PartAtOnce",
    "PartAtOnce",
    "PartAtOnce",
    "PartAtOnce",
    "PartAtOnce",
    "PartAtOnce",
]
color_table = [
    "",
    "red",
    "green",
    "orange",
    "green",
    "purple",
    "purple",
    "purple",
    "purple",
    "purple",
    "purple",
]
fig, ax = plt.subplots()

df = pd.read_csv(name, sep=";")
if var == "d":
    toplot = df.drop(
        df[
            (df["n"] != int(sys.argv[3]))
            | (df["t"] != int(sys.argv[4]))
            | (df["ss_type"] != sys.argv[5])
        ].index
    )
elif var == "t":
    toplot = df.drop(
        df[
            (df["n"] != int(sys.argv[3]))
            | (df["d"] != int(sys.argv[4]))
            | (df["ss_type"] != sys.argv[5])
        ].index
    )
for key, grp in toplot.groupby("algo"):
    temp = grp.groupby([var], as_index=False).agg(
        min=pd.NamedAgg(column="approx_part", aggfunc="min"),
        max=pd.NamedAgg(column="approx_part", aggfunc="max"),
    )
    grp.groupby([var]).mean("approx_part").plot(
        ax=ax, y="approx_part", label=label_table[int(key)], color=color_table[int(key)]
    )

    ax.fill_between(
        x=var,
        y1="min",
        y2="max",
        data=temp,
        color=mpl.colors.to_rgba(color_table[int(key)], 0.05),
    )

temp = toplot.groupby([var], as_index=False).agg(
    min=pd.NamedAgg(column="approx_random", aggfunc="min"),
    max=pd.NamedAgg(column="approx_random", aggfunc="max"),
    approx_random=pd.NamedAgg(column="approx_random", aggfunc="mean"),
)
temp.groupby([var]).mean("approx_random").plot(
    ax=ax, y="approx_random", label="random", color="blue"
)

ax.fill_between(
    x=var, y1="min", y2="max", data=temp, color=mpl.colors.to_rgba("blue", 0.05)
)

# plt.plot(tlist, [1/m.sqrt(t) for t in tlist], c='grey', linestyle='dotted', label="$\\frac{1}{\\sqrt{t}}$")
# plt.plot(tlist, [1/m.sqrt(t)**((int(sys.argv[2])+1)/(int(sys.argv[2]))) for t in tlist], c='grey', linestyle='dashdot', label="$\\frac{1}{\\sqrt{t^{\\frac{d+1}{d}}}}$")
if var == "d":
    plt.plot(
        dlist,
        [1 / m.sqrt(int(sys.argv[4])) for d in dlist],
        c="grey",
        linestyle="dotted",
        label="$\\frac{1}{\\sqrt{t}}$",
    )
    plt.plot(
        dlist,
        [1 / m.sqrt(int(sys.argv[4])) ** ((d + 1) / d) for d in dlist],
        c="grey",
        linestyle="dashdot",
        label="$\\frac{1}{\\sqrt{t^{\\frac{d+1}{d}}}}$",
    )
elif var == "t":
    plt.plot(
        tlist,
        [1 / m.sqrt(t) for t in tlist],
        c="grey",
        linestyle="dotted",
        label="$\\frac{1}{\\sqrt{t}}$",
    )
    plt.plot(
        tlist,
        [1 / m.sqrt(t) ** ((2 + 1) / 2) for t in tlist],
        c="grey",
        linestyle="dashdot",
        label="$\\frac{1}{\\sqrt{t^{\\frac{d+1}{d}}}}$",
    )

ax.set_xlabel(var, fontsize=25)
ax.set_ylabel("epsilon", fontsize=25)

plt.legend(fontsize=20)
plt.show()

