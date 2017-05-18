from general import Satellite, Vector
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def createOutput():
    creation_time = '42'
    global ResultList
    ResultList = open('OrbitOutput_' + creation_time + '.txt', 'w')
    # ResultList.write('x\ty\tz\n')


def orbitOutput(t, r: Vector, V: Vector):
    ResultList.write(str(r.x) + ' ' + str(r.y) + ' ' + str(r.z) + '\n')
    print(t)


def finalize():
    global ResultList
    ResultList.close()

    ResultList = open('OrbitOutput_42.txt', 'r')

    x = []
    y = []
    z = []
    ResultList.readline()
    for line in ResultList:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))
        z.append(float(line.split()[2]))

    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    # plt.plot(x, y)
    # plt.show()

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(x, y, z)
    plt.show()

    ResultList.close()



