import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def is_valid_interval(a, b, solution):
    solution_in_interval = False
    if sp.Interval(a, b).contains(solution):
        solution_in_interval = True
        
    if not solution_in_interval:
        return False
    return True

def get_nodes(a, b, m):
    nodes = np.linspace(a, b, m)
    return nodes

def compute_omega(x, nodes, k):
    n = len(nodes) - 1

    omega = 1
    for j in range (n + 1):
        if j == k:
            continue
        omega *= (x - nodes[j])
    
    return omega

def lagrange_interpolation(x, nodes):
    n = len(nodes) - 1
    
    L = 0
    for k in range(n + 1):
        numerator, denominator = 1, 1
        for j in range (n + 1):
            if j == k:
                continue
            
            numerator *= (x - nodes[j])
            denominator *= (nodes[k] - nodes[j])
        
        L += f_func(nodes[k]) * (numerator / denominator)

    return L

def forward_interpolation(nodes):
    #if not sp.is_monotonic(f, sp.Interval.open(a, b)):
        x = sp.Symbol('x')
    
        L = 0
        for k in range(n + 1):
            numerator, denominator = 1, 1
            for j in range (n + 1):
                if j == k:
                    continue
                
                numerator *= (x - nodes[j])
                denominator *= (nodes[k] - nodes[j])
            
            L += f_func(nodes[k]) * (numerator / denominator)


        root = sp.nsolve(L, x, 0)
        return root, L

def backward_intepolation(nodes):
    #if not sp.is_monotonic(f, sp.Interval.open(a, b)):
        y = sp.Symbol('y')
    
        L = 0
        for k in range(n + 1):
            numerator, denominator = 1, 1
            for j in range (n + 1):
                if j == k:
                    continue
                
                numerator *= (y - f_func(nodes[j]))
                denominator *= (f_func(nodes[k]) - f_func(nodes[j]))
            
            L += nodes[k] * (numerator / denominator)
    
        L_func = sp.lambdify(y, L)
        root = L_func(0) 
        return root, L

# Пошук точного розв'язку
x = sp.Symbol('x')
f = np.e ** x - 2 * (x - 1) ** 2
#f = sp.ln(x) + x - 2
f_func = sp.lambdify(x, f)
print(f"Функція: {f}")

exact_solution = sp.nsolve(f, x, 1)
print(f"Точний розв'язок: {exact_solution}")

# Дослідження інтервалу
a = -1
b = 1

if not is_valid_interval(a, b, exact_solution):
   print("Жоден розв'язок не належить вказаному інтервалу")
   exit()

m = 10 # к-сть вузлів інтерполяції
n = m - 1 # степінь полінома

# Отримання рівновіддалених вузлів інтерполяції
nodes = get_nodes(a, b, m)
print(f"Вузли інтерполяції: {nodes}\nКрок інтерполяції: {nodes[1] - nodes[0]}")

# Пряма інтерполяція
forward_solution, forward_polinom = forward_interpolation(nodes)
print(f"Пряма інтерполяція: {forward_solution}")

# Зворотня інтерполяція
backward_solution, backward_polinom = backward_intepolation(nodes)
print(f"Зворотня інтерполяція: {backward_solution}")

# Графічне зображення функції та поліномів
x_vals = np.linspace(-4, 4, 100)
y_vals = f_func(x_vals)

# Обмеження розмірів
plt.figure(figsize=(10, 5))
plt.xlim(-2, 2)
plt.ylim(-10, 10)

# Вісі координат
plt.axhline(0, color='black', lw=1)
plt.axvline(0, color='black', lw=1)

# Графік функції
plt.plot(x_vals, y_vals, label='f(x)', color='blue')

# Графік прямої інтерполяції
forward_polinom_func = sp.lambdify(x, forward_polinom)
forward_polinom_vals = forward_polinom_func(x_vals)
plt.plot(x_vals, forward_polinom_vals, label='Пряма інтерполяція', color='orange', linestyle='--')

# Графік зворотної інтерполяції
backward_polinom_func = sp.lambdify(sp.Symbol('y'), backward_polinom)
backward_polinom_vals = backward_polinom_func(x_vals)
plt.plot(x_vals, backward_polinom_vals, label='Зворотня інтерполяція', color='green', linestyle='--')

# Точки інтерполяції
plt.scatter(nodes, f_func(nodes), color='black', label='Вузли інтерполяції')

# Розв'язки
plt.scatter(forward_solution, 0, color='red', label='Корінь прямої інтерполяції')
plt.scatter(backward_solution, 0, color='purple', label='Корінь зворотної інтерполяції')
plt.scatter(exact_solution, 0, color='brown', label='Точний розв\'язок')

plt.legend()
plt.grid()
plt.show()