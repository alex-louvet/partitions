import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np
import math as m

df = pd.read_csv("experiments_grid.csv", sep=";")
toplot = df.drop(df[(df["n"] != int(sys.argv[1])) | (df["d"] != int(sys.argv[2])) |(df["t"] != int(sys.argv[3])) | (df["ss_type"] != sys.argv[4])].index)
#toplot = df.drop(df[(df["n"] != int(sys.argv[1])) | (df["p"] != float(sys.argv[2])) | (df["t"] != float(sys.argv[3]))].index)

''' exp on k
print(toplot)
print("\multicolumn{1}{c|}{"+sys.argv[1]+","+sys.argv[2]+","+sys.argv[3]+"} & " + str(toplot.iat[0,7]) + "& \multicolumn{1}{c|}{" + str('{:.4}'.format(toplot.iat[0,11])) + "} & " + str(toplot.iat[1,7]) + " & \multicolumn{1}{c|}{" + str('{:.4}'.format(toplot.iat[1,11])) + "}&" + str(toplot.iat[2,7]) + " & \multicolumn{1}{c}{" + str('{:.4}'.format(toplot.iat[2,11])) + "}\\\\")
'''

#print(toplot)

cross = toplot.groupby('algo').mean("max_crossing")
rt = toplot.groupby('algo').mean("runtime")
print(cross)
print(rt)

''' additional exp rate c
for i in range(5):
    print("\multicolumn{1}{c|}{" + str(0.5*i) + "} & " + str('{:.5}'.format(cross.iat[5+i,5])) + " & \multicolumn{1}{c|}{" + str('{:.5}'.format(cross.iat[5+i,8])) + "} & " + str('{:.5}'.format(cross.iat[i,5])) + " & \multicolumn{1}{c|}{" + str('{:.5}'.format(cross.iat[i,8])) + "} & " + str('{:.5}'.format(cross.iat[10+i,5])) + " & \multicolumn{1}{c|}{" + str('{:.5}'.format(cross.iat[10+i,8])) + "}")
'''

''' rate violation
print("\multicolumn{1}{c|}{"+sys.argv[1]+","+sys.argv[2]+","+sys.argv[3]+"} & \multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[1,8])) + "} & \multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[0,8])) + "}& \multicolumn{1}{c}{" + str('{:.4}'.format(cross.iat[2,8])) + "}\\\\")
'''
'''
print("\multicolumn{1}{c|}{"+sys.argv[1]+","+sys.argv[2]+","+sys.argv[3]+"} & \multicolumn{1}{c|}{" + str('{:.4}'.format(2*int(sys.argv[3])**(1-1/(int(sys.argv[2]))))) + "} &" + str('{:.4}'.format(cross.iat[1,5])) + "& \multicolumn{1}{c|}{" + str('{:.4}'.format(rt.iat[1,9])) + "} & " + str('{:.4}'.format(cross.iat[0,5])) + " & \multicolumn{1}{c|}{" + str('{:.4}'.format(rt.iat[0,9])) + "}&" + str('{:.4}'.format(cross.iat[2,5])) + " & \multicolumn{1}{c|}{" + str('{:.4}'.format(rt.iat[2,9])) + "} & " + str('{:.4}'.format(cross.iat[3,5])) + " & \multicolumn{1}{c}{" + str('{:.4}'.format(rt.iat[3,9])) + "}\\\\")
'''

#power law
'''
print("\multicolumn{1}{c|}{"+sys.argv[1]+","+sys.argv[2]+","+sys.argv[3]+"} & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} &"  + str('{:.4}'.format(cross.iat[0,5])) + " & \multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[0,9])) + "}&" + str('{:.4}'.format(cross.iat[2,5])) + " & \multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[2,9])) + "} & " + str('{:.4}'.format(cross.iat[3,5])) + " & \multicolumn{1}{c}{" + str('{:.4}'.format(cross.iat[3,9])) + "}\\\\")
#print("\multicolumn{1}{c|}{"+sys.argv[1]+","+sys.argv[2]+","+sys.argv[3]+"} & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} &"  + str('{:.4}'.format(cross.iat[0,5])) + " & \multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[0,9])) + "}&" + str('{:.4}'.format(cross.iat[2,5])) + " & \multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[2,9])) + "} & " + str('{:.4}'.format(cross.iat[3,5])) + " & \multicolumn{1}{c}{" + str('{:.4}'.format(cross.iat[3,9])) + "}\\\\")
'''
'''osm
print("\multicolumn{1}{c|}{"+sys.argv[2]+"," + sys.argv[3] + "} &" + str(toplot.iat[1,7]) + "& \multicolumn{1}{c|}{" + str('{:.6}'.format(toplot.iat[1,11])) + "} & " + str(toplot.iat[0,7]) + " & \multicolumn{1}{c|}{" + str('{:.5}'.format(toplot.iat[0,11])) + "}&" + str(toplot.iat[2,7]) + " & \multicolumn{1}{c|}{" + str('{:.5}'.format(toplot.iat[2,11])) + "} & " + str(toplot.iat[3,7]) + " & \multicolumn{1}{c}{" + str('{:.4}'.format(toplot.iat[3,11])) + "}\\\\")
'''

#print("\multicolumn{1}{c|}{"+sys.argv[2]+"," + sys.argv[3] + "} &" + str(cross.iat[1,5]) + "& \multicolumn{1}{c|}{" + str('{:.6}'.format(cross.iat[1,9])) + "} & " + str(cross.iat[0,5]) + " & \multicolumn{1}{c|}{" + str('{:.5}'.format(cross.iat[0,9])) + "}&" + str(cross.iat[2,5]) + " & \multicolumn{1}{c}{" + str('{:.5}'.format(cross.iat[2,9])) + "}\\\\")


print("\multicolumn{1}{c|}{"+sys.argv[1]+","+sys.argv[2]+","+sys.argv[3]+"} & \multicolumn{1}{c|}{" + str('{:.4}'.format(2*int(sys.argv[3])**(1-1/(int(sys.argv[2]))))) + "} &" + str('{:.4}'.format(cross.iat[0,5])) + "&" + str('{:.4}'.format(cross.iat[0,6])) + "& \multicolumn{1}{c|}{" + str('{:.4}'.format(rt.iat[0,9])) + "} & " + str('{:.4}'.format(cross.iat[1,5])) + "&" + str('{:.4}'.format(cross.iat[1,6])) + " & \multicolumn{1}{c|}{" + str('{:.4}'.format(rt.iat[1,9])) + "}&" + str('{:.4}'.format(cross.iat[3,5])) + "&" + str('{:.4}'.format(cross.iat[3,6])) + " & \multicolumn{1}{c|}{" + str('{:.4}'.format(rt.iat[3,9])) + "} & " + str('{:.4}'.format(cross.iat[4,5])) + "&" + str('{:.4}'.format(cross.iat[4,6])) + " & \multicolumn{1}{c}{" + str('{:.4}'.format(rt.iat[4,9])) + "}\\\\")

''' approx
print("\multicolumn{1}{c|}{" + sys.argv[1]+"," +sys.argv[2]+"," + sys.argv[3] + "} & \multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[0,11]/3 + cross.iat[1,11]/3 + cross.iat[2,11]/3)) + "} & \multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[0,10])) + "} & \multicolumn{1}{c|}{" + str('{:.4}'.format(cross.iat[1,10])) + "} & \multicolumn{1}{c}{" + str('{:.4}'.format(cross.iat[2,10])) + "}\\\\")
'''