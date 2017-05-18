from general import Vector
from math import sqrt, sin, cos, factorial, atan2
from numpy import zeros
from scipy import special


def load_geo_coefficients():
    data = open("C_and_S_data.txt", "r")
    global C
    global S
    global J
    C = zeros([51, 51])
    S = zeros([51, 51])
    J = zeros([51])
    for line in data:
        n = int(line.split()[0])
        m = int(line.split()[1])
        if n == 51:
            break
        Cv = float(line.split()[2])
        Sv = float(line.split()[3])
        C[n, m] = Cv * sqrt(factorial(n - m) * (2 * n + 1) * (2 - int(m == 0)) / factorial(n + m))
        S[n, m] = Sv * sqrt(factorial(n - m) * (2 * n + 1) * (2 - int(m == 0)) / factorial(n + m))
        if m == 0:
            J[n] = Cv * sqrt(factorial(n - m) * (2 * n + 1) * (2 - int(m == 0)) / factorial(n + m))
    data.close()


def forceConstructor(r:Vector, V:Vector, time):
    if time == 0:
        load_geo_coefficients()

    forceSummary = Vector(0, 0, 0)
    forceSummary += acceleration_main(r) + geopotential_J(r, time, 20) + geopotential_C_S(r, time, 20)

    return forceSummary


def acceleration_main(r: Vector):
    Mu = 398600.4418e9
    a = - Mu * r / r.ro**3
    return a


# Функція простого поліному Лежандра заданого порядку
def Pn(n, value):
    return special.lpmv(0, n, value)


# Функція для отримання значення першої похідної поліному Лежандра
def P_prime_n(n, value):
    return (1 + n) / (1 - value**2) * (- Pn(n + 1, value) + value * Pn(n, value))


# Функція для отримання значення приєднаного поліному Лежандра
def Pnm(n, m, value):
    return special.lpmv(m, n, value)


# Функція для отримання значення похідної приєднаного поліному Лежандра множиного на косинус
def cosPnm_prime(n, m, sinFi, cosFi):
    return - n * sinFi / cosFi * Pnm(n, m, sinFi) + (n + m) / cosFi * Pnm(n-1, m, sinFi)


def geopotential_J(r: Vector, time, N=2):
    Re = 6378.13630e3
    Mu = 398600.4418 * 1e9
    theta = 7292115e-11 * time

    x =   r.x * cos(theta) + r.y * sin(theta)
    y = - r.x * sin(theta) + r.y * cos(theta)
    z =   r.z

    ro = sqrt(x ** 2 + y ** 2 + z ** 2)
    sinFi = float(z / ro)
    cosFi = float(sqrt(x ** 2 + y ** 2) / ro)
    sinLambda = float(y / sqrt(x ** 2 + y ** 2))
    cosLambda = float(x / sqrt(x ** 2 + y ** 2))

    aJr = 0
    aJfi = 0

    for n in range(N):
        aJr += J[n] * (Re / ro) ** n * (n + 1) * Pn(n, sinFi)
        aJfi -= J[n] * (Re / ro) ** n * cosFi * P_prime_n(n, sinFi)

    aJr *= - Mu / ro ** 2
    aJfi *= - Mu / ro ** 2

    axe = - sinFi * cosLambda * aJfi + cosFi * cosLambda * aJr
    aye = - sinFi * sinLambda * aJfi + cosFi * sinLambda * aJr
    aze = cosFi * aJfi + sinFi * aJr

    ax = axe * cos(theta) - aye * sin(theta)
    ay = aye * cos(theta) + axe * sin(theta)
    az = aze

    return Vector(ax, ay, az)


def geopotential_C_S(r: Vector, time, N=2):
    Re = 6378.13630e3
    Mu = 398600.4418 * 1e9

    theta = 7292115e-11 * time

    x = r.x * cos(theta) + r.y * sin(theta)
    y = - r.x * sin(theta) + r.y * cos(theta)
    z = r.z

    ro = sqrt(x**2 + y**2 + z**2)
    sinFi = float(z / ro)
    cosFi = float(sqrt(x**2 + y**2) / ro)
    sinLambda = float(y / sqrt(x**2 + y**2))
    cosLambda = float(x / sqrt(x**2 + y**2))
    Lambda = atan2(y, x)

    aCSr = 0
    aCSlambda = 0
    aCSfi = 0

    for n in range(N):
        for m in range(1, n + 1):
            aCSr -= (Re / ro)**n * (n + 1) * Pnm(n, m, sinFi) * \
                    (C[n, m] * cos(m * Lambda) + S[n, m] * sin(m * Lambda))

            aCSlambda += (Re / ro)**n * m / cosFi * Pnm(n, m, sinFi) * \
                         (- C[n, m] * sin(m * Lambda) + S[n, m] * cos(m * Lambda))

            aCSfi += cosPnm_prime(n, m, sinFi, cosFi) * (C[n, m] * cos(m * Lambda) + S[n, m] * sin(m * Lambda))

    aCSr *= - Mu / ro**2
    aCSlambda *= - Mu / ro**2
    aCSfi *= - Mu / ro**2

    axe = cosFi * cosLambda * aCSr - sinLambda * aCSlambda - sinFi * cosLambda * aCSfi
    aye = cosFi * sinLambda * aCSr + cosLambda * aCSlambda - sinFi * sinLambda * aCSfi
    aze = sinFi * aCSr + cosFi * aCSfi

    ax = axe * cos(theta) - aye * sin(theta)
    ay = aye * cos(theta) + axe * sin(theta)
    az = aze

    return Vector(ax, ay, az)