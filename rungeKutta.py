from general import Satellite, Vector
from forces import forceConstructor
from math import sqrt
import outputMananger

class Integrator:
    def __init__(self, step, steptype='static'):
        self.info = '''Runge Kutta integrator'''
        self.h = step
        self.h_new = step
        # steptype may be static or dynamic. Dynamic requires a one additional summation action,
        # but the step size won't be equally separated in time
        self.steptype = steptype

        if self.steptype != 'static' and self.steptype != 'dynamic':
            print('''ERROR. Ordertype may be only 'static' or 'dynamic'. Simulation continues with static type.''')
            self.steptype = 'static'

        # Here one sets coefficients for RK. They may be found anywhere in the literature.

        # self.c = [0, 0, 1/2, 1/2, 1]
        # self.a = [[0],
        #           [0],
        #           [0, 1/2],
        #           [0, 0, 1/2],
        #           [0, 0, 0, 1]]
        # self.b = [0, 1/6, 1/3, 1/3, 1/6]
        # self.b_star = [0, 1/6, 1/3, 1/3, 1/6]

        self.c = [0, 0, 1 / 10, 2 / 9, 3 / 7, 3 / 5, 4 / 5, 1, 1]
        self.a = [[0],
                  [0],
                  [0, 1 / 10],
                  [0, -2 / 81, 20 / 81],
                  [0, 615 / 1372, -270 / 343, 1053 / 1372],
                  [0, 3243 / 5500, -54 / 55, 50949 / 71500, 4998 / 17875],
                  [0, -26492 / 37125, 72 / 55, 2808 / 23375, -24206 / 37125, 338 / 459],
                  [0, 5561 / 2376, -35 / 11, -24117 / 31603, 899983 / 200772, -5225 / 1836, 3925 / 4056],
                  [0, 465467 / 266112, -2945 / 1232, -5610201 / 14158144, 10513573 / 3212352, -424325 / 205632,
                   376225 / 454272, 0]]
        self.b = [0, 61 / 864, 0, 98415 / 321776, 16807 / 146016, 1375 / 7344, 1375 / 5408, -37 / 1120, 1 / 10]
        self.b_star = [0, 821/10800, 0, 19683/71825, 175273/912600, 395/3672, 785/2704, 3/50, 0]

        # *star parameters are for dynamic mode only. if they are not set -- the simulation will be static

        if not self.b_star:
            print('''b_star == None. Simulation continues with static type.''')
            self.steptype = 'static'

        if self.steptype == 'dynamic':
            self.epsilon = 1e-10
            self.rtol = self.epsilon # relative error
            self.atol = self.epsilon # absolute error
            self.err = 1
            self.S = 0.9

    def integrate(self, s: Satellite, maxTime):
        t = 0
        outputMananger.createOutput()
        while t <= maxTime:
            # compute values of r and V on the next step
            s.r, s.V = self.integrateStep(t, s.r, s.V)
            # use ``proxy'' to output the data about satellite time, position, etc.
            outputMananger.orbitOutput(t, s.r, s.V)
            t += self.h
        # finalize proxy work (plot orbit, close files, etc.)
        outputMananger.finalize()

    def integrateStep(self, time, r: Vector, V: Vector):
        # Just RK classical routine.

        k = [Vector(0, 0, 0)] * len(self.b)
        l = [Vector(0, 0, 0)] * len(self.b)
        self.h = self.h_new

        for i in range(1, len(self.b)):
            k[i] = self.h * (V + sum(self.a[i][j] * l[j] for j in range(1, i)))
            l[i] = self.h * forceConstructor(r + sum(self.a[i][j] * k[j] for j in range(1, i)),
                                             V + sum(self.a[i][j] * l[j] for j in range(1, i)),
                                             time)

        if self.steptype == 'dynamic':
            r_star = r + sum(self.b_star[i] * k[i] for i in range(1, len(self.b_star)))
            V_star = V + sum(self.b_star[i] * l[i] for i in range(1, len(self.b_star)))

        r += sum(self.b[i] * k[i] for i in range(1, len(self.b)))
        V += sum(self.b[i] * l[i] for i in range(1, len(self.b)))

        if self.steptype == 'dynamic':
            delta = [r.x - r_star.x, V.x - V_star.x,
                     r.y - r_star.y, V.y - V_star.y,
                     r.z - r_star.z, V.z - V_star.z]
            scale = [self.atol + abs(r.x) * self.rtol, self.atol + abs(r.x) * self.rtol,
                     self.atol + abs(r.y) * self.rtol, self.atol + abs(r.y) * self.rtol,
                     self.atol + abs(r.z) * self.rtol, self.atol + abs(r.z) * self.rtol]
            err = sqrt(1 / 2 * sum((delta[i] / scale[i])**2 for i in range(6)))
            self.h_new = self.h * self.S * (1 / err)**(1 / 5)
            # print('------------>', self.h_new)

        return r, V
