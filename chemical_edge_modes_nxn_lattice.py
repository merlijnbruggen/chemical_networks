import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig
import random
# nxn lattice with 4 internal sites
def construct_markov_chain_on_lattice(n):
    """
    Returns the adjacency matrix for an nxn lattice as descirbed in Tang et al (2021)
    """

    #####
    # Internal transitions
    #####
    
    P_internal = np.zeros((4*(n**2),4*(n**2)))

    for j in range(0, n):
        basis_vector = []

        for i in range(0, 2*n):
            basis_vector.append(i)
        basis_vector = [entry + j*4*n for entry in basis_vector]
        

        for idx in basis_vector[::2]:
            
            P_internal[idx, idx+1] = 1
            P_internal[idx+1, idx+1+len(basis_vector)] = 1
            P_internal[idx+1+len(basis_vector), idx+len(basis_vector)] = 1
            P_internal[idx+len(basis_vector), idx] = 1

    
    plt.imshow(P_internal)
    plt.title("P_internal")
    plt.show()

    #####
    # External transitions
    #####
    P_external = np.zeros((4*(n**2),4*(n**2)))
    

    #VERTICAL external transitions

    # initialize basis vector for vertical external transitions
    # the range is until n-1 because the bottom row does not have any vertical external transitions
    
    basis_vector = []

    # initialize basis vector for vertical transitions
    
    for j in range(0, n-1):
        basis_vector = []
        for _ in range(2*n, 4*n):
            basis_vector.append(_)
        
        basis_vector = [i+j*4*n for i in basis_vector]
        

        for idx in basis_vector[1::2]:
            P_external[idx, idx+2*n] = 1
        for idx in basis_vector[::2]:
            P_external[idx+2*n, idx] = 1
        
        

    #######
    #HORIZONTAL external transitions
    #######

    # initialize basis vector for horizontal tranisitions
    
    for j in range(0, n-1):
        basis_vector = []
        for _ in range(1, 4*n**2 - n + 1, 2*n):
            basis_vector.append(_)  #initialize horizontal basis vector
        
        basis_vector = [entry + j*2 for entry in basis_vector]
        


        for idx in basis_vector[1::2]:
            P_external[idx + 1, idx] = 1
            
        for idx in basis_vector[::2]:
            P_external[idx, idx + 1] = 1
        
            

    plt.imshow(P_external)
    plt.title("P_external")
    plt.show()
   

    return P_internal, P_external


n = 2
P_internal, P_external = construct_markov_chain_on_lattice(n)


gamma_int = 1
gamma_ext = 100

P_total = gamma_int * P_internal + gamma_ext * P_external

# We add a small diagonal matrix to the transition matrix
# I expect this to be physically relevant bc it is conceivable that sometimes no transition happens (just as in a gillespie algorithm right?)
P_total = P_total + 0 * np.identity(P_total.shape[0])

initial_state=random.randint(0, 4*(n**2))


'''
def construct_transition_matrix(S):

    # Get the number of nodes (rows/columns of the stoichiometric matrix)
    n = S.shape[0]
    
    # Initialize the transition matrix with zeros
    P = np.zeros((n, n))
    
    # Iterate over each row of the stoichiometric matrix
    for i in range(n):
        # Calculate the sum of absolute flux magnitudes out of node i (this will normalize the row)
        # This to account for the difference in the internal and external transition probabilities
        flux_sum = np.sum(np.abs(S[i, :]))
        
        # If there's no flux, skip to the next row (transition probability remains 0)
        if flux_sum == 0:
            print(f"there is not flux out of node {i}")
            continue
        
        # Normalize each element of the row by the total outgoing flux
        for j in range(n):
            P[i, j] = np.abs(S[i, j]) / flux_sum

    # After constructing the matrix, we ensure that each row sums to 1
    # Check and normalize rows where sum != 1
    for row in range(n):
        row_sum = np.sum(P[row, :])  # Sum the elements in the current row
        
        # If the row sum isn't 1, adjust the diagonal element to ensure the row sums to 1
        if row_sum != 1:
            # If there is flux, adjust the diagonal to make the row sum 1
            P[row, row] += 1 - row_sum

    
    
    plt.imshow(P)
    plt.title("transition matrix for the Markov chain")
    plt.show()
    

    
    return P
transition_matrix = construct_transition_matrix(P_total)



def simulate_markov_chain(P, initial_state, n_steps):
    print(f"intitial state = {initial_state}")
    n_states = P.shape[0]  # Number of states
    states = [initial_state]

    for _ in range(n_steps):
        current_state = states[-1]
        next_state = np.random.choice(n_states, p=P[current_state])
        states.append(next_state)
    print(f"simlutated {states=}")
    return states

n_steps = 500
trajectory = simulate_markov_chain(transition_matrix, initial_state=initial_state, n_steps=500)


import matplotlib.animation as animation
def visualize_markov_chain(trajectory, n_rows, n_cols):
    # Set up the figure and axis
    fig, ax = plt.subplots(figsize=(6, 6))
    
    # Create a scatter plot for the grid points
    grid_x, grid_y = np.meshgrid(np.arange(4*n_cols), np.arange(4*n_rows))
    ax.scatter(grid_y, n_rows - 1 - grid_x, s=10, c='b', label="Grid Points")  # Plot grid points
    
    # Create a scatter plot to represent the current state on the grid
    scatter = ax.scatter([], [], s=50, c='r', label="Current State")
    
    # Set up grid limits and labels
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ax.set_xticks(np.arange(0, n_cols, 1))
    ax.set_yticks(np.arange(0, n_rows, 1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    # Animation update function
    def update(frame):
        x, y = divmod(trajectory[frame], n_cols)
        scatter.set_offsets([y, n_rows - 1 - x])
        return scatter,
    
    # Create the animation
    ani = animation.FuncAnimation(fig, update, frames=len(trajectory), interval=200, blit=True)
    
    
    # Add grid point and animation legend
    ax.legend(loc="upper left")
    plt.title(f"Markov Chain on Grid. {gamma_int=} and {gamma_ext=}. Grid size: {n} by {n}")
    plt.show()
    
    #save animation
    #writervideo = animation.FFMpegWriter(fps=2) 
    #ani.save('chemical_network_edge_modes.mp4', writer=writervideo) 
    plt.close() 

# Visualize the running of the Markov chain on the grid (6x8 grid)
visualize_markov_chain(trajectory, n_rows=n*2, n_cols=n*2)


def compute_occupancy_distribution(trajectory, lattice_size, initial_state):
    """
    Takes a trajectory of a Markov chain on a lattice and compute where the systems spends time
    """

    counter = np.zeros(lattice_size**2)
    
    for idx in trajectory:
        counter[idx] += 1
    
    plt.bar(np.arange(lattice_size**2), counter, edgecolor='black')
    plt.xticks(np.arange(lattice_size**2), labels=np.arange(lattice_size**2))
    plt.title(f"distribution of state occupancy. {n_steps=}, {gamma_int=}, {gamma_ext=}")
    plt.ylabel("counts")
    plt.xlabel("states in Markov chain")
    #plt.show()
    
    plt.clf()
    reshaped_counter = counter.reshape((lattice_size,lattice_size))
    plt.imshow(reshaped_counter)
    plt.colorbar()
    plt.title(f"heatmap of occupancy for a {n}x{n} lattice. {n_steps=}, {gamma_int=}, {gamma_ext=}, {initial_state=}")
    #plt.savefig(f"{gamma_ext}.png")
    plt.show()
    
    return counter

compute_occupancy_distribution(trajectory, lattice_size=2*n, initial_state=initial_state)



def compute_stationary_distribution(P, lattice_size):
    eigenvals, eigenvecs = np.linalg.eig(P.T)

    idx = np.isclose(eigenvals, 1)

    stationary_dist = abs(eigenvecs[:, idx])
    

    stationary_dist_reshaped = stationary_dist.reshape((lattice_size, lattice_size))
    
    plt.imshow(stationary_dist_reshaped)
    
    plt.title(f"stationary distribution for a {n}x{n} lattice. {gamma_int=}, {gamma_ext=}")
    plt.colorbar()
    plt.show()
    return stationary_dist

compute_stationary_distribution(transition_matrix, lattice_size=2*n)
'''



'''
for _ in np.arange(1, 15, 0.5):
    gamma_ext = _ 
    P_total = gamma_int * P_internal + gamma_ext * P_external
    transition_matrix = construct_transition_matrix(P_total)
    n_steps = 200
    initial_state=36
    trajectory = simulate_markov_chain(transition_matrix, initial_state=36, n_steps=200)
    compute_occupancy_distribution(trajectory, lattice_size=2*n, initial_state=initial_state)
'''


######
# Use a master equation to learn more about the stationary distribution of the system
# The solutions of the ODEs all diverge, it doesn't seem very workable
######

from scipy.integrate import odeint
from scipy.integrate import solve_ivp

def solve_system_of_odes(n, P_internal, P_external, P_total):
    # Step 1: Compute the degree matrix
    adjacency_matrix = P_internal + P_external
    degree_matrix = np.zeros((4 * (n**2), 4 * (n**2)))
    
    for i in range(adjacency_matrix.shape[0]):
        degree = np.sum(adjacency_matrix[i, :])
        degree_matrix[i, i] = degree
        

    # Step 2: Compute the Laplacian matrix W
    W = P_total - degree_matrix
    print(f"{W=}")
    # Step 3: Define the system of ODEs (the master equation)
    def dpdt(p, t, W):
        print(f"{np.dot(W,p)=}")
        return np.dot(W, p)

    # Step 4: Set initial conditions (random in this case)
    initial_conditions = np.zeros(4 * (n**2))
    for i in range(len(initial_conditions) - 1):
        initial_conditions[i] = random.uniform(1, 2)

    # Step 5: Define the time grid
    t = np.linspace(0, 0.5, 1000)  # Adjust the time range as needed

    # Step 6: Solve the system using odeint
    solution = odeint(dpdt, initial_conditions, t, args=(W,))

    # Step 7: Plot the solution
    plt.figure(figsize=(10, 6))

    # Loop through each function (column of the solution matrix)
    for i in range(solution.shape[1]):
        plt.plot(t, solution[:, i], label=f"Node {i+1}")

    # Add labels and title
    plt.xlabel('Time')
    plt.ylabel('Concentration of node')
    plt.title('Solutions of the Master Equation')

    # Add legend
    plt.legend()

    # Show the plot
    plt.show()

#solve_system_of_odes(n, P_internal, P_external, P_total)


#####
# Now try to simulate the lattice using a Gillespie algorithm
#####


# Define concentrations for each node
concentrations = [[0] for _ in range(4 * (n**2))]
concentrations[0][0] = 1

'''
non_zero_internal = [np.nonzero(P_internal[:, i])[0] for i in range(P_internal.shape[1])]
non_zero_external = [np.nonzero(P_external[:, i])[0] for i in range(P_external.shape[1])]

def rates(P_internal, P_external, concentrations, non_zero_internal, non_zero_external):
    rates = np.zeros(P_internal.shape[0])  # Use numpy for faster array operations
    num_nodes = P_internal.shape[0]
    
    for i in range(num_nodes):
        # INTERNAL
        if len(non_zero_internal[i]) > 0:
            non_zero_idx = non_zero_internal[i][0]
            rates[i] += P_internal[non_zero_idx, i] * concentrations[non_zero_idx][-1]
        
        # EXTERNAL
        if len(non_zero_external[i]) > 0:
            non_zero_idx = non_zero_external[i][0]
            rates[i] += P_external[non_zero_idx, i] * concentrations[non_zero_idx][-1]
    
    for j in range(num_nodes):
        # INTERNAL
        if len(non_zero_internal[j]) > 0:
            non_zero_idx = non_zero_internal[j][0]
            rates[j] += -P_internal[j, non_zero_idx] * concentrations[non_zero_idx][-1]
        
        # EXTERNAL
        if len(non_zero_external[j]) > 0:
            non_zero_idx = non_zero_external[j][0]
            rates[j] += -P_external[j, non_zero_idx] * concentrations[non_zero_idx][-1]
    
    return rates
'''

# the rates array is computed incorrectly. For each node, there should be TWO events: depletion and addition. Now there is only one entry per node. 

def rates(P_internal, P_external, concentrations):
    rates_gain = [0] * P_internal.shape[0]
    rates_loss = [0] * P_internal.shape[0]
    num_nodes = P_internal.shape[0]
    
    
    # Rows of P_internal and P_external
    # popoulate rates_gain
    for i in range(num_nodes):
        # INTERNAL: non-zero column indices
        column = P_internal[:, i]
        non_zero_indices = np.nonzero(column)[0]    # there is only 1 non zero entry (no backwards reactions)
        if non_zero_indices.size > 0:
            non_zero_idx = non_zero_indices[0]
            rates_gain[i] += gamma_int * concentrations[non_zero_idx][-1]
            

        # EXTERNAL: non-zero column indices
        column_ext = P_external[:, i]
        non_zero_indices = np.nonzero(column_ext)[0]
        if non_zero_indices.size > 0:
            non_zero_idx = non_zero_indices[0]
            rates_gain[i] += gamma_ext * concentrations[non_zero_idx][-1]
            

    
    # Columns of P_internal and P_external
    # The events associated with these are the removal of a unit of concentration
    # populate rates_loss
    for j in range(num_nodes):
        # INTERNAL: non-zero row indices
        row = P_internal[j, :]
        non_zero_indices = np.nonzero(row)[0]
        if non_zero_indices.size > 0:
            non_zero_idx = non_zero_indices[0]
            rates_loss[j] += gamma_int * concentrations[non_zero_idx][-1]
          
        # EXTERNAL: non-zero row indices
        row_ext = P_external[j, :]
        non_zero_indices = np.nonzero(row_ext)[0]
        if non_zero_indices.size > 0:
            non_zero_idx = non_zero_indices[0]
            rates_loss[j] += gamma_ext * concentrations[non_zero_idx][-1]
            
    # Interleave the arrays
    rates_total = np.empty(len(rates_gain) + len(rates_loss))
    rates_total[0::2] = rates_gain  # Fill even indices 
    rates_total[1::2] = rates_loss  # Fill odd indices

    print(f"{rates_total=}, {len(rates_total)=}")
    return rates_total



time = [0]
t_final = 1

while time[-1] < t_final:
    
    '''
    # get the current rates
    rate_values = rates(P_internal, P_external, concentrations)
    print(f"{rate_values=}")
    # When will the next event happen?
    rate_sum = sum(rate_values)
    print(f"{rate_sum=}")
    tau = np.random.exponential(scale=1/rate_sum)
    time.append(time[-1] + tau)

    # Which event will happen?
    # List all the possible events
    rand = random.uniform(0,1)

    
    if rand * rate_sum < rate_values[0]:
        concentrations[0] = np.concatenate([concentrations[0], [concentrations[0][-1] + 1]])

        # the rest doesn't change
        for idx in range(1, len(concentrations) - 1):
            concentrations[idx] = np.concatenate([concentrations[idx], concentrations[idx]])

    elif rand * rate_sum > rate_values[0] and rand * rate_sum < sum(rate_values[:2]):
        concentrations[0] = np.concatenate([concentrations[0], [concentrations[0][-1] - 1]])

        for idx in range(1, len(concentrations) - 1):
            concentrations[idx] = np.concatenate([concentrations[idx], concentrations[idx]])
    '''

    
    # Compute rates
    rate_values = rates(P_internal, P_external, concentrations)
    
    # When will the next event happen?
    rate_sum = sum(rate_values)
    tau = 1/rate_sum * np.log(1/(np.random.random()))
    time.append(time[-1] + tau)

    # Select the next event
    rand = random.uniform(0, 1)
    cumulative_rates = np.cumsum(rate_values)
    print(f"{cumulative_rates=}")
    selected_event = np.searchsorted(cumulative_rates, rand * rate_sum)
    print(f"{selected_event=}")
    #TODO: somehow the concentrations become negative at longer times and the selected event is always event 8. I don't know why
    
    
    
    node = selected_event % P_internal.shape[0] # the lattice size
    loss_or_gain = selected_event % 2   # gain (mod2 =0) or loss (mod2 =1)
    print(f"{node=}, {loss_or_gain=}")

    if loss_or_gain == 0:
        concentrations[node][-1] += 1
    elif loss_or_gain == 1:
        print(f"LOSING CONCENTRATION")
        concentrations[node][-1] -= 1
        

    #and all the other nodes stay the same
    for i in range(len(concentrations)):
        if i != selected_event % P_internal.shape[0]:
            concentrations[i].append(concentrations[i][-1])
            print(f"{i=} and {concentrations[i][-1]=}")

    '''
    # Update concentrations for the selected event
    # TODO: is the event selection correct?
    if rand * rate_sum < cumulative_rates[selected_event % P_internal.shape[0]] and selected_event % 2 == 0:
        concentrations[selected_event].append(concentrations[selected_event][-1] + 1)
    elif rand * rate_sum < cumulative_rates[selected_event % P_internal.shape[0]] and selected_event % 2 == 1:
        concentrations[selected_event].append(concentrations[selected_event][-1] - 1)
        print("concentration decreases!")
    '''
    '''
    elif concentrations[selected_event][-1] > cumulative_rates[selected_event] and concentrations[selected_event][-1] < cumulative_rates[selected_event + 1] and :
        # TODO: This elif statement never gets executed 
        concentrations[selected_event].append(concentrations[selected_event][-1] - 1)
        '''
    '''
    # Update concentrations for unchanged nodes
    for i in range(len(concentrations)):
        if i != selected_event % P_internal.shape[0]:
            concentrations[i].append(concentrations[i][-1])
    '''
    
    

    
    print(f"time: {time[-1]}")
    print(f"{rate_sum=}")



    





    