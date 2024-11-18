import numpy as np
import sympy as sp

def chebishev_zeros(a, b):
    zeros = []
    for k in range (n):
        zeros.append((a + b) / 2 + (b - a) / 2 * np.cos((2 * k + 1) * np.pi / (2 * n))) # нулі полінома Чебишева
    return zeros

x = sp.Symbol('x')
f = np.e ** x - 2 * (x - 1) ** 2
exact_solution = sp.nsolve(f, x, 1)
print(f"Exact solution: {exact_solution}")

m = 10 # к-сть вузлів інтерполяції
n = m - 1 # степінь полінома