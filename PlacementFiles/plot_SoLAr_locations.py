import matplotlib.pyplot as plt
import csv
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


## Read in the SiPM locations (We want to plot them in green)
with open('./SiPM', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    x = []
    y = []
    z = []
    color = []
    for row in reader:
        x.append(float(row[1]))
        y.append(float(row[2]))
        z.append(float(row[3]))
        color.append('g')

with open('./Pixel', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    for row in reader:
        x.append(float(row[1]))
        y.append(float(row[2]))
        z.append(float(row[3]))
        color.append('r')


fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)


length = int(len(x))

for i in range(length):
    if( color[i] == 'r'):
        # Add the squares at the pixel locations (red)
        # Note that this xoors, ycoords is only for the z=const plane.
        # Other planes = other orientations (see 1_....py)
        xsize = 0.3/2
        ysize = 0.3/2
        zsize = 0.0
        x_low = x[i] - xsize
        x_high = x[i] + xsize
        y_low = y[i] - ysize
        y_high = y[i] + ysize
        xcoords  =  [ x_low, x_low, x_high, x_high ]
        ycoords  =  [ y_low, y_high, y_high, y_low ]
        zcoords  =  [ z[i], z[i], z[i], z[i] ]
    if( color[i] == 'g'):
        # Add the squares at the SiPM locations (green)
        xsize = 0.6/2
        ysize = 0.6/2
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

    xcoords = [0 , 0, 1, 1]
    ycoords = [0 , 1, 1, 0]
    zcoords = [0 , 0, 0, 0]


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
# ax.set_xlim(0, 300)
# ax.set_ylim(0, 300)
# ax.set_zlim(0, 300)

ax.auto_scale_xyz([0, 50], [0, 50], [0, 1])


plt.show()
