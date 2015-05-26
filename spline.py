import numpy
import sympy

class SplineInterpolate:
    def __init__(self, x_list, y_list, boundary_condition='second', const=(0, 0)):
        n = len(x_list)
        xy_list = list(zip(x_list, y_list))
        xy_list.sort()
        self.xy = xy_list
        dx = []
        for i in range(n - 1):
            dx.append(xy_list[i + 1][0] - xy_list[i][0])
        mu = [0, ]
        lam = [0, ]
        for i in range(1, n - 1):
            mu.append(dx[i - 1] / (dx[i - 1] + dx[i]))
            lam.append(dx[i] / (dx[i - 1] + dx[i]))
        deviation = []
        for i in range(1, n - 1):
            deviation.append(6 * (get_deviation(xy_list, i) - get_deviation(xy_list, i - 1)) / (dx[i - 1] + dx[i]))
        a = numpy.zeros((n, n))
        if boundary_condition == 'second':
            mu.append(0)
            deviation.insert(0, 2 * const[0])
            deviation.append(2 * const[1])
        else:
            lam[0] = 1
            mu.append(1)
            deviation.insert(0, 6 / dx[0] * (get_deviation(xy_list, 0) - const[0]))
            deviation.append(6 / dx[-1] * (const[1] - get_deviation(xy_list, -2)))
        for i in range(a.shape[0]):
            a[i][i] = 2
            a[i][i - 1] = mu[i] if i != 0 else 0
            try:
                a[i][i + 1] = lam[i]
            except IndexError:
                pass
        m = numpy.linalg.solve(a, deviation)
        m = list(m)
        self.m = m
        self.dx = dx
        self._get_coeffs()

    def get_knots(self):
        return list(zip(*self.xy))[0]

    def _get_coeffs(self):
        self.coeffs = []
        x = sympy.symbols('x')
        for i in range(len(self.xy) - 1):
            polynomial = self.m[i] * (self.xy[i + 1][0] - x) ** 3 / (6 * self.dx[i])
            polynomial += self.m[i + 1] * (x - self.xy[i][0]) ** 3 / (6 * self.dx[i])
            polynomial += (self.xy[i][1] - self.m[i] * self.dx[i] ** 2 / 6) * (self.xy[i + 1][0] - x) / self.dx[i]
            polynomial += (self.xy[i + 1][1] - self.m[i + 1] * self.dx[i] ** 2 / 6) * (x - self.xy[i][0]) / self.dx[i]
            if polynomial == 0:
                self.coeffs.append([0, 0, 0, 0])
            else:
                self.coeffs.append(sympy.poly(polynomial).all_coeffs())

    def get_coeffs(self):
        return self.coeffs

    def __call__(self, *args, **kwargs):
        x = args[0]
        if hasattr(x, '__contains__'):
            ans = []
            for x1 in x:
                ans.append(self(x1))
            return ans
        for i in range(1, len(self.xy) - 1):
            if x < self.xy[i][0]:
                return get_polynomial_value(self.coeffs[i-1], x)
        return get_polynomial_value(self.coeffs[-1], x)


def get_deviation(xy_list, i):
    xy1 = xy_list[i + 1]
    xy0 = xy_list[i]
    return (xy1[1] - xy0[1]) / (xy1[0] - xy0[0])

def get_polynomial_value(coeffs, x):
    y = 0
    for i, coeff in enumerate(coeffs[::-1]):
        y += coeff * x ** i
    return y
