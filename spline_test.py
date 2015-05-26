from spline import SplineInterpolate

X_LIST = [0.25, 0.30, 0.39, 0.45, 0.53]
Y_LIST = [0.5000, 0.5477, 0.6245, 0.6708, 0.7280]

f = SplineInterpolate(X_LIST, Y_LIST, boundary_condition='first', const=(1, 0.6868))
print(f.get_knots())
print(f.get_coeffs())
print(f(0.25))

f = SplineInterpolate(X_LIST, Y_LIST, boundary_condition='second', const=(0, 0))
print(f.get_knots())
print(f.get_coeffs())
print(f(0.25))
