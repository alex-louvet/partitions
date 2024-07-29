import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

name = sys.argv[1]

name2 = sys.argv[2]


def listify(x):
    x = x[:-1].split(",")
    temp = []
    for y in x:
        temp.append(int(y))
    return np.array(temp)

df = pd.read_csv(name, sep=";")
df['table'] = df['table'].apply(listify)

toplot = df.drop(df[(df["n"] != int(sys.argv[1])) | (df["d"] != int(sys.argv[2])) |(df["t"] != int(sys.argv[3])) | (df["algo"] != 1)].index)

sum = np.array([0]*(int(sys.argv[3])))
for index,row in toplot.iterrows():
    sum += np.array(row['table'])
temp = []
for x in sum:
    temp.append(x/100)
sum = np.array(temp)

toplot = df.drop(df[(df["n"] != int(sys.argv[1])) | (df["d"] != int(sys.argv[2])) |(df["t"] != int(sys.argv[3])) | (df["algo"] != 4)].index)

sum2 = np.array([0]*(int(sys.argv[3])))
for index,row in toplot.iterrows():
    sum2 += np.array(row['table'])
temp = []
for x in sum2:
    temp.append(x/100)
sum2 = np.array(temp)

df2 = pd.read_csv(name2, sep=";")
df2['table'] = df2['table'].apply(listify)

toplot = df2.drop(df2[(df2["n"] != int(sys.argv[1])) | (df2["d"] != int(sys.argv[2])) |(df2["t"] != int(sys.argv[3])) | (df2["algo"] != 1)].index)

sum3 = np.array([0]*(int(sys.argv[3])))
for index,row in toplot.iterrows():
    sum3 += np.array(row['table'])
temp = []
for x in sum3:
    temp.append(x/100)
sum3 = np.array(temp)

toplot = df2.drop(df2[(df2["n"] != int(sys.argv[1])) | (df2["d"] != int(sys.argv[2])) |(df2["t"] != int(sys.argv[3])) | (df2["algo"] != 4)].index)

sum4 = np.array([0]*(int(sys.argv[3])))
for index,row in toplot.iterrows():
    sum4 += np.array(row['table'])
temp = []
for x in sum4:
    temp.append(x/100)
sum4 = np.array(temp)

print(max(sum2))

#plt.pcolor([sum,sum2], cmap="viridis", edgecolor="w", vmin=0, vmax=max(np.max(sum),np.max(sum2)), linewidths=2)
c = plt.pcolor([sum,sum2,sum3,sum4], cmap="Greens", edgecolor="w", vmin=0, vmax=max(1,max(sum2)), linewidths=2)
plt.colorbar(c,shrink=0.1, pad=0.01,extend='neither')
plt.axis('equal')
plt.axis('off')
fs = 10
plt.text(-1,0,"4",fontsize=fs)
plt.text(-1,1,"3",fontsize=fs)
plt.text(-1,2,"2",fontsize=fs)
plt.text(-1,3,"1",fontsize=fs)
plt.show()