import matplotlib.pyplot as plt
import math as m

file1 = open("example.txt", "r")
var = False
for l in file1.readlines():
    if not var:
        pts = l.split(";")[:-1]
        for i in range(len(pts)):
            pts[i] = pts[i].split(",")
            if len(pts[i]) == 3:
                for j in range(len(pts[i])):
                    pts[i][j] = float(pts[i][j])
    else:
        root = l.split(",")
        for j in range(len(root)):
                root[j] = float(root[j])
    
    
    if var:
        plt.scatter([x[0] for x in pts if x[2] != 0],[x[1] for x in pts if x[2] != 0],c=[x[2] for x in pts if x[2] != 0], cmap='gray')
        plt.scatter([root[0]],[root[1]], color='red')
        plt.show()
        plt.clf()
    var = not var