from math import sqrt, radians, cos, sin

class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.ro = sqrt(x**2 + y**2 + z**2)

    def __add__(self, other):
        if other == 0:
            return self
        else:
            return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __mul__(self, other):
        if type(self) == type(other):
            return self.x * other.x + self.y * other.y + self.z * other.z
        else:
            return Vector(self.x * other, self.y * other, self.z * other)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return Vector(self.x / other, self.y / other, self.z / other)

    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)

    def __neg__(self):
        return Vector(- self.x, - self.y, - self.z)

    def __str__(self):
        return str(self.x) + '\t' + str(self.y) + '\t' + str(self.z)

    def cross(self, other):
        return Vector(self.y * other.z - self.z * other.y,
                      self.z * other.x - self.x * other.z,
                      self.x * other.y - self.y * other.x)

    def __abs__(self):
        return Vector(abs(self.x), abs(self.y), abs(self.z))

    def max(self):
        if self.x >= self.y and self.x >= self.z:
            return self.x

        if self.y >= self.x and self.y >= self.z:
            return self.y

        if self.z >= self.y and self.z >= self.x:
            return self.z



class Satellite:
    def __init__(self, m, x, y, z, Vx, Vy, Vz, JulianCenchuries=0, system='car'):
        self.m = m
        if system == 'ecl':
            epsilon = radians(23.4392916666667 - 0.013004167 * JulianCenchuries - 1.66667E-07 * JulianCenchuries**2 +
                              5.02778E-07 * JulianCenchuries**3)
            self.r = Vector(x, y * cos(epsilon) - z * sin(epsilon), z * cos(epsilon) + y * sin(epsilon))
            self.V = Vector(Vx, Vy * cos(epsilon) - Vz * sin(epsilon), Vz * cos(epsilon) + Vy * sin(epsilon))
        else:
            self.r = Vector(x, y, z)
            self.V = Vector(Vx, Vy, Vz)

    def __str__(self):
        return 'r = ' + str(self.r * 1e-3) + ' km\n' + \
               'V = ' + str(self.V * 1e-3) + ' km/s\n'