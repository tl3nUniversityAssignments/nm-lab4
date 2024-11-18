import numpy as np
import sympy as sp

m = 10 # к-сть вузлів інтерполяції
n = 3 # степінь полінома

def chebishev_zeros(a, b):
    zeros = []
    for k in range (n):
        zeros.append((a + b) / 2 + (b - a) / 2 * np.cos((2 * k + 1) * np.pi / (2 * n))) # нулі полінома Чебишева
    return zeros

print(chebishev_zeros(-1, 1))