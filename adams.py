from general import Satellite, Vector
from forces import forceConstructor
import outputMananger
from adams_coef_generator import create_coeff
import os.path
from numpy import copy

class Integrator:
    def __init__(self, step, ordertype='static', order=5):
        self.info = '''Adams integrator'''
        self.h = step

        # Checking if coefficients are available
        if not os.path.isfile('Adams_Predictor.txt') or not os.path.isfile('Adams_Corrector.txt'):
            create_coeff()
        self.coeff_list = get_coefficients(40, 'Adams_Predictor.txt')
        self.coeff_list_star = get_coefficients(41, 'Adams_Corrector.txt')

        self.currentOrder = 1
        self.order = order
        self.ordertype = ordertype

        # ordertype may be static or dynamic. Dynamic requires 6 more computations of the right side
        if self.ordertype != 'static' and self.ordertype != 'dynamic':
            print('''ERROR. Ordertype may be only 'static' or 'dynamic'. Simulation continues with static type.''')
            self.ordertype = 'static'

        self.V_list = []
        self.a_list = []

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
        self.a_list.append(forceConstructor(r, V, time))
        self.V_list.append(V)

        if self.ordertype == 'static':
            alpha = copy(self.coeff_list[self.currentOrder])
            alpha_star = copy(self.coeff_list[self.currentOrder + 1])

            pV = V + sum(alpha[i] * self.a_list[i] for i in range(len(alpha))) * self.h
            pr = r + sum(alpha[i] * self.V_list[i] for i in range(len(alpha))) * self.h
            self.V_list.insert(0, pV)
            self.a_list.insert(0, forceConstructor(pr, pV, time))

            V += sum(alpha_star[i] * self.a_list[i] for i in range(len(alpha_star))) * self.h
            r += sum(alpha_star[i] * self.V_list[i] for i in range(len(alpha_star))) * self.h

            self.V_list[0] = V
            self.a_list[0] = forceConstructor(r, V, time)

            if self.currentOrder < self.order:
                self.currentOrder += 1

            if len(self.a_list) > 50:
                self.a_list.pop()
                self.V_list.pop()

        if self.ordertype == 'dynamic':
            r_cu, V_cu = self.estimate2(time, r, V, self.currentOrder, self.currentOrder + 1)
            r_ce, V_ce = self.estimate2(time, r, V, self.currentOrder, self.currentOrder)
            sigma_c = abs(r_cu - r_ce).max()

            r_mi, V_mi = self.estimate2(time, r, V, self.currentOrder, self.currentOrder)
            r_me, V_me = self.estimate2(time, r, V, self.currentOrder, self.currentOrder - 1)
            sigma_m = abs(r_mi - r_me).max()

            r_pl, V_pl = self.estimate2(time, r, V, self.currentOrder + 1, self.currentOrder + 2)
            r_pe, V_pe = self.estimate2(time, r, V, self.currentOrder + 1, self.currentOrder + 1)
            sigma_p = abs(r_pl - r_pe).max()

            r, V = r_cu, V_cu
            if sigma_c > sigma_p and self.currentOrder < 30:
                self.currentOrder += 1
                r, V = r_pl, V_pl
            elif sigma_c > sigma_m and self.currentOrder > 1:
                self.currentOrder -= 1
                r, V = r_mi, V_mi
            print(self.currentOrder)

            self.V_list.insert(0, V)
            self.a_list.insert(0, forceConstructor(r, V, time))

            if len(self.a_list) > 50:
                self.a_list.pop()
                self.V_list.pop()



        return r, V

    def estimate2(self, time, r, V, p_order, c_order):
        alpha = copy(self.coeff_list[p_order])
        alpha_star = copy(self.coeff_list_star[c_order])
        a_list = self.a_list
        V_list = self.V_list

        pV = V + sum(alpha[i] * self.a_list[i] for i in range(len(alpha))) * self.h
        pr = r + sum(alpha[i] * self.V_list[i] for i in range(len(alpha))) * self.h
        V_list.insert(0, pV)
        a_list.insert(0, forceConstructor(pr, pV, time))

        V += sum(alpha_star[i] * self.a_list[i] for i in range(len(alpha_star))) * self.h
        r += sum(alpha_star[i] * self.V_list[i] for i in range(len(alpha_star))) * self.h

        return r, V

def get_coefficients(k, path):
    k += 1
    data = open(path, 'r')
    ksi = [[0.0]]
    j = 0
    for line in data:
        ksi.append([])
        j += 1
        r = line[2:].split()
        for i in r:
            ksi[j].append(float(i))
        s = line.split()
        if int(s[0]) == k:
            break
    return ksi