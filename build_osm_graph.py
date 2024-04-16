import sys
import json
import matplotlib.pyplot as plt

#json osm file
file1 = open("data/roads-paris.json", "r")
#name of file to export the edges
file2 = open("paris-edges.txt","w")
# name of file to export the coordinates of each node
file3 = open("paris-geo.txt","w")
data = json.load(file1)

trans = {}
c = 1
edges = []
for x in data['elements']:
    ide = str(x['geometry'][0]['lat'])+","+str(x['geometry'][0]['lon'])
    ide2 = str(x['geometry'][-1]['lat'])+","+str(x['geometry'][-1]['lon'])
    if not ide in trans and not ide2 in trans:
        trans[ide] = c
        file3.write(str(x['geometry'][0]['lat']*1000)+" "+str(x['geometry'][0]['lon']*1000) + '\n')
        c += 1
        trans[ide2] = c
        file3.write(str(x['geometry'][-1]['lat']*1000)+" "+str(x['geometry'][-1]['lon']*1000) + '\n')
        c += 1
        edges.append([ide2, ide])
    elif not ide in trans and ide2 in trans:
        trans[ide] = c
        file3.write(str(x['geometry'][0]['lat']*1000)+" "+str(x['geometry'][0]['lon']*1000) + '\n')
        c += 1
        edges.append([ide2, ide])
    elif not ide2 in trans and ide in trans:
        trans[ide2] = c
        file3.write(str(x['geometry'][-1]['lat']*1000)+" "+str(x['geometry'][-1]['lon']*1000) + '\n')
        c += 1
        edges.append([ide2, ide])
    else:
        edges.append([ide2, ide])
    

for e in edges:
    ''' drawing
    s0 = e[0].split(",")
    s1 = e[1].split(",")
    a = [float(s0[0]),float(s0[1])]
    b = [float(s1[0]),float(s1[1])]
    plt.scatter([a[0],b[0]],[a[1],b[1]],color='red')
    plt.plot([a[0],b[0]],[a[1],b[1]],color='black')
    '''

    file2.write(str(trans[e[0]]) + " " + str(trans[e[1]]) + '\n')

plt.show()
