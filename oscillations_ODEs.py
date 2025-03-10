import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def hill(H, h, x):
    print(f"hill={(H**h) / (H**h + x**h)}")
    return (H**h) / (H**h + x**h + 1e-8)

gamma_int = 1
gamma_ext = 1
h = 1   #Hill steepness
H = 1
h_ext = 2
total_mass = []
epsilon = 1
def func(t, x):
    """returns the system of ODEs that correspond
    to our model."""
    
    # First node (node A)
    dx0dt = gamma_int * x[3] * x[0] / (1 + x[0])  - gamma_int * x[1] * x[0] / (1 + x[0]) + gamma_ext * x[6] / (1 + x[0]) 
    dx1dt = gamma_int * x[0] * x[1] / (1 + x[1]) - gamma_int * x[2] * x[1] / (1 + x[1])
    dx2dt = gamma_int * x[1] * x[2] / (1 + x[2]) - gamma_int * x[3] * x[2] / (1 + x[2]) - gamma_ext * x[2] / (1 + x[4]) 
    dx3dt = gamma_int * x[2] * x[3] / (1 + x[3]) - gamma_int * x[0] * x[3] / (1 + x[3])
   
    # Second node (node B)
    dx4dt = gamma_int * x[7] * x[4] / (1 + x[4]) - gamma_int * x[5] * x[4] / (1 + x[4]) + gamma_ext * x[2] / (1 + x[4]) 
    dx5dt = gamma_int * x[4] * x[5] / (1 + x[5]) - gamma_int * x[6] * x[5] / (1 + x[5])
    dx6dt = gamma_int * x[5] * x[6] / (1 + x[6]) - gamma_int * x[7] * x[6] / (1 + x[6]) - gamma_ext * x[6] / (1 + x[0])
    dx7dt = gamma_int * x[6] * x[7] / (1 + x[7]) - gamma_int * x[4] * x[7] / (1 + x[7])
    

    total_mass = (np.sum([x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]]))
    
    
    print(f"{t=}")
    print(f"{total_mass=}")

    if np.max([dx0dt, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt]) > 100:
        print(f"derivative greater than 100 encountered: {np.max([dx0dt, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt])}")
    return [dx0dt, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt]

t_span = (0, 50)
t = np.linspace(0, 50, 10000)
y0 = [1,epsilon,epsilon,epsilon,1,epsilon,epsilon,epsilon]

#y0 = np.random.uniform(0.1, 1.5, size=8)
#y0 = np.ones(8)
solution = solve_ivp(func, t_span, y0=y0, t_eval=t)
print(f"Final time reached: {solution.t[-1]}")






plt.plot(solution.t, solution.y[0], label='x0')
plt.plot(solution.t, solution.y[1], label='x1')
plt.plot(solution.t, solution.y[2], label='x2')
plt.plot(solution.t, solution.y[3], label='x3')
plt.plot(solution.t, solution.y[4], label='x4')
plt.plot(solution.t, solution.y[5], label='x5')
plt.plot(solution.t, solution.y[6], label='x6')
plt.plot(solution.t, solution.y[7], label='x7')
plt.xlabel("Time (s)")
plt.ylabel("solutions")
plt.title("Solutions of chemical network")
plt.legend()
plt.show()

