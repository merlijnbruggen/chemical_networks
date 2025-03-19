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
epsilon = 0.0001
def func(t, x):
    """returns the system of ODEs that correspond
    to our model."""
    
    '''
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
    '''

    #2nd try
    dx0dt = -gamma_int * x[0] + gamma_int * x[3] + gamma_ext * x[0] * x[6] #+ gamma_ext * x[6] / (1 + x[0]) #
    dx1dt = -gamma_int * x[1] + gamma_int * x[0]
    dx2dt = -gamma_int * x[2] + gamma_int * x[1] - gamma_ext * x[2] * x[4] #- gamma_ext * x[2] / (1 + x[4]) #
    dx3dt = -gamma_int * x[3] + gamma_int * x[2]

    dx4dt = -gamma_int * x[4] + gamma_int * x[7] + gamma_ext *x[4] * x[2]##+ gamma_ext * x[2] / (1 + x[4]) #
    dx5dt = -gamma_int * x[5] + gamma_int * x[4]
    dx6dt = -gamma_int * x[6] + gamma_int * x[5] - gamma_ext * x[6] * x[0]# #- gamma_ext * x[6] / (1 + x[0]) # 
    dx7dt = -gamma_int * x[7] + gamma_int * x[6]


    return [dx0dt, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt]

t_span = (0, 25)
t = np.linspace(0, 25, 100000)
#y0 = [2,2, gamma_int * 2/(gamma_int + gamma_ext * 2), gamma_int * 2/(gamma_int + gamma_ext * 2), 2, 2, gamma_int * 1/(gamma_int + gamma_ext * 1),gamma_int * 1/(gamma_int + gamma_ext * 1)]
#y0 = [0, 0, 0, 0, 2, 2, 2, 2]  #fixed point
#y0 = [-1+epsilon,-1+epsilon,-1+epsilon,-1+epsilon, 0, 0, 0, 0] #fixed point
y0 = np.random.uniform(0.1, 1.5, size=8)
#y0 = np.ones(8)
solution = solve_ivp(func, t_span, y0=y0, t_eval=t)
print(f"Final time reached: {solution.t[-1]}")





from adjustText import adjust_text
plt.figure(figsize=(8, 6))
texts = []  # Store text objects

for i in range(8):
    plt.plot(solution.t, solution.y[i], label=f'x{i}')
    text = plt.text(solution.t[-1], solution.y[i][-1], f'x{i}', fontsize=9)
    texts.append(text)  # Collect text objects

# Automatically adjust text positions to avoid overlap
adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray'))
plt.xlabel("Time (s)")
plt.ylabel("solutions")
plt.title("Solutions of 2 cell model")
plt.legend()
plt.show()

plt.plot(solution.y[1], solution.y[4])
plt.title(f"phase space of x[1] vs. x[4]")
plt.show()

total_mass = np.sum(solution.y, axis=0)
total_mass_nodeA = np.sum(solution.y[0:3], axis=0)
total_mass_nodeB = np.sum(solution.y[4:7], axis=0)
plt.plot(solution.t, total_mass)
plt.plot(solution.t, total_mass_nodeA)
plt.plot(solution.t, total_mass_nodeB)
plt.xlabel('Time')
plt.ylabel('Total mass per node')
plt.title('Total mass of each node and total mass')
plt.legend()
plt.show()