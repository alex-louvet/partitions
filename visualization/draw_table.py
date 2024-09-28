import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np
import math as m

df = pd.read_csv("results.csv", sep=";")
if sys.argv[5] != "pl":
    toplot = df.drop(df[(df["n"] != int(sys.argv[1])) | (df["d"] != int(sys.argv[2])) |(df["t"] != int(sys.argv[3])) | (df["ss_type"] != sys.argv[4])].index)
    cross = toplot.groupby('algo').mean("max_crossing")
    rt = toplot.groupby('algo').mean("runtime")
else:
    toplotpl = df.drop(df[(df["n"] != int(sys.argv[1])) | (df["p"] != float(sys.argv[2])) | (df["t"] != float(sys.argv[3]))| (df["ss_type"] != sys.argv[4])].index)
    crosspl = toplotpl.groupby('algo').mean("max_crossing")

#rate violation
if sys.argv[5] == "rate":
    print("\\multicolumn{1}{c|}{"+sys.argv[1]+","+sys.argv[2]+","+sys.argv[3]+"} & \\multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[1,8])) + "} & \\multicolumn{1}{c}{" + str('{:.4}'.format(cross.iat[0,8])) + "}\\\\")

vctable = [5.8,4,3,6.5,4.4,3]

#power law
if sys.argv[5] == "pl":
    print("\\multicolumn{1}{c|}{"+sys.argv[1]+","+sys.argv[2]+","+sys.argv[3]+"} & \\multicolumn{1}{c|}{" + str(vctable[(int(sys.argv[1]) > 5000)*3+int(float(sys.argv[2])*2)-4]) + "} & \\multicolumn{1}{c|}{" + str('{:.3}'.format(int(sys.argv[3])**(1-1/(vctable[(int(sys.argv[1]) > 5000)*3+int(float(sys.argv[2])*2)-4])))) + "} &"  + str('{:.4}'.format(crosspl.iat[0,5])) + " & \\multicolumn{1}{c|}{" + str('{:.4}'.format(crosspl.iat[0,9])) + "}&" + str('{:.4}'.format(crosspl.iat[1,5])) + " & \\multicolumn{1}{c}{" + str('{:.4}'.format(crosspl.iat[1,9])) + "}\\\\")

#grid
if sys.argv[5] == "grid":
    print("\\multicolumn{1}{c|}{"+sys.argv[1]+","+sys.argv[2]+","+sys.argv[3]+"} & \\multicolumn{1}{c|}{" + str('{:.4}'.format(2*int(sys.argv[3])**(1-1/(int(sys.argv[2]))))) + "} &" + str('{:.4}'.format(cross.iat[0,5])) + "&" + str('{:.4}'.format(cross.iat[0,6])) + "& \\multicolumn{1}{c|}{" + str('{:.4}'.format(rt.iat[0,9])) + "} & " + str('{:.4}'.format(cross.iat[1,5])) + "&" + str('{:.4}'.format(cross.iat[1,6])) + " & \\multicolumn{1}{c|}{" + str('{:.4}'.format(rt.iat[1,9])) + "}&" + str('{:.4}'.format(cross.iat[2,5])) + "&" + str('{:.4}'.format(cross.iat[2,6])) + " & \\multicolumn{1}{c}{" + str('{:.4}'.format(rt.iat[2,9])) + "}\\\\")

#approx
if sys.argv[5] == "approx":
    print("\\multicolumn{1}{c|}{" + sys.argv[1]+"," +sys.argv[2]+"," + sys.argv[3] + "} & \\multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[0,11]/3 + cross.iat[1,11]/3 + cross.iat[2,11]/3)) + "} & \\multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[0,10])) + "} & \\multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[1,10])) + "} & \\multicolumn{1}{c}{" + str('{:.4}'.format(cross.iat[2,10])) + "}\\\\")

#facebook
if sys.argv[5] == "fb":
    print("\\multicolumn{1}{c|}{" +sys.argv[2]+"," + sys.argv[3] + "} & " + str('{:.4}'.format(cross.iat[0,5])) + " & \\multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[0,9])) + "} & " + str('{:.4}'.format(cross.iat[1,5])) + " & \\multicolumn{1}{c}{" + str('{:.4}'.format(cross.iat[1,9])) + "}\\\\")