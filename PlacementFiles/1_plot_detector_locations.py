from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
import csv

with open('./full_dune.txt', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    x = []
    y = []
    z = []
    orientation = []
    n = 0
    for row in reader:
        x.append(float(row[1]))
        y.append(float(row[2]))
        z.append(float(row[3]))
        if(float(row[5]) == 1):
            orientation.append('r')
        elif(float(row[5]) == 2):
            orientation.append('b')
        elif(float(row[5]) == 3):
            orientation.append('g')
        else:
            print("Error: orientation not found")

# with open('./SoLAr_1cm.txt', 'r') as csvfile:
#     reader = csv.reader(csvfile, delimiter=' ')
#     for row in reader:
#         x.append(float(row[1]))
#         y.append(float(row[2]))
#         z.append(float(row[3]))
#         if(float(row[5]) == 1):
#             orientation.append('r')
#         elif(float(row[5]) == 2):
#             orientation.append('b')
#         elif(float(row[5]) == 3):
#             orientation.append('g')
#         else:
#             print("Error: orientation not found")



from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


fig = plt.figure()

ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)


length = int(len(x))
print(length)

for i in range(length):

    if( orientation[i] == 'r'):
        xsize = 0.0
        ysize = 3
        zsize = 3
        y_low = y[i] - ysize
        y_high = y[i] + ysize
        z_low = z[i] - zsize
        z_high = z[i] + zsize
        xcoords  =  [ x[i], x[i], x[i], x[i]]
        ycoords  =  [ y_low,y_high,y_high,y_low]
        zcoords  =  [ z_low,z_low,z_high,z_high]

    if( orientation[i] == 'b'):
        xsize = 3
        ysize = 0.0
        zsize = 3
        x_low = x[i] - xsize
        x_high = x[i] + xsize
        z_low = z[i] - zsize
        z_high = z[i] + zsize
        xcoords  =  [ x_low, x_low, x_high, x_high ]
        ycoords  =  [ y[i], y[i], y[i], y[i] ]
        zcoords  =  [ z_low, z_high, z_high, z_low ]

    if( orientation[i] == 'g'):
        xsize = 2.
        ysize = 2.
        zsize = 0.0
        x_low = x[i] - xsize
        x_high = x[i] + xsize
        y_low = y[i] - ysize
        y_high = y[i] + ysize
        xcoords  =  [ x_low, x_low, x_high, x_high ]
        ycoords  =  [ y_low, y_high, y_high, y_low ]
        zcoords  =  [ z[i], z[i], z[i], z[i] ]

    verts = [list(zip(xcoords,ycoords,zcoords))]
    if(orientation[i] == 'g'):
        ax.add_collection3d(Poly3DCollection(verts, facecolor=orientation[i], linewidths=1, edgecolors='k', alpha=.25))
    else:
        ax.add_collection3d(Poly3DCollection(verts, facecolor=orientation[i], linewidths=1, edgecolors='k', alpha=.25))

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
# ax.set_xlim(0, 300)
# ax.set_ylim(0, 300)
# ax.set_zlim(0, 300)

ax.auto_scale_xyz([0, 100], [0, 100], [0, 100])


plt.show()
