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

    '''
    plt.imshow(P_internal)
    plt.title("P_internal")
    plt.show()
    '''
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
        
            
    '''
    plt.imshow(P_external)
    plt.title("P_external")
    plt.show()
    '''
   

    return P_internal, P_external


n = 4
P_internal, P_external = construct_markov_chain_on_lattice(n)

gamma_int = 1
gamma_ext = 1

P_total = gamma_int * P_internal + gamma_ext * P_external

plt.imshow(P_total)
plt.title("P_total")
plt.show()


# I am not even using this function in my gillespie loop
def rates(P_total, node):
    """
    Given a node, return the rate of exiting that node. Since there are no reverse reactions, one always
    has to leave a given node
    """

    rates = []
    for entry in P_total[:,node]:
        rates.append(entry)
    rates.sort()

    return rates

num_nodes = num_nodes = P_total.shape[0]
time = [0]
tmax = 50
initial_condition = np.random.choice(num_nodes)
print(f"{initial_condition=}")
node = initial_condition
trajectory = []

while time[-1] < tmax:
    rate_values = []
    non_zero_indices = np.nonzero(P_total[node, :])
    
    for i in non_zero_indices:
        rate_values.append(P_total[node, i])

    

    # time to next event
    rate_sum = np.sum(rate_values)
    tau = np.random.exponential(scale=1/rate_sum) 
    time.append(time[-1] + tau)
    
    # next event
    rand = random.uniform(0,1)
    cumulative_rates = np.cumsum(P_total[node, :])
    

    selected_event = np.searchsorted(cumulative_rates, rand * rate_sum)

    node = selected_event
    trajectory.append(node)
    
    print(f"{time[-1]=}")




import matplotlib.animation as animation
def visualize_reaction(trajectory, n_rows, n_cols, dynamic_playback):
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
        # Update the title with the current "time" value
        ax.set_title(f"Gillespie simulation on grid. {n_rows} x {n_cols}. Time: {time[frame]:.2f}", fontsize=14)
        return scatter, ax

    # Create the animation
    ani = animation.FuncAnimation(fig, update, frames=len(trajectory), interval=100, blit=False)
    # Add grid point and animation legend
    ax.legend(loc="upper left")
    #plt.title(f"Gillespie simulation on grid. {n_rows} x {n_cols}.")
    plt.show()
    

    if dynamic_playback == True:
        # Have the frame update be equal to the gillespie random times
        # Calculate frame intervals based on the time differences
        time_intervals = np.diff(time) * 1000 # Convert to milliseconds
        time_intervals = np.append(time_intervals, time_intervals[-1])  # Extend the last frame's interval
            # Function to dynamically control the interval
        def frame_interval_gen():
            for interval in time_intervals:
                yield interval
        
        # Create the animation
        ani = animation.FuncAnimation(
            fig, 
            update, 
            frames=len(trajectory), 
            interval=1,  # Placeholder, dynamically set below
            blit=True
        )
        
        # Update the frame intervals dynamically
        ani.event_source.stop()  # Stop initial playback
        interval_gen = frame_interval_gen()
        def update_interval(*args):
            try:
                ani.event_source.interval = next(interval_gen)
            except StopIteration:
                ani.event_source.stop()
        
        ani.event_source.add_callback(update_interval)
        ani.event_source.start()  # Start playback with dynamic intervals
        
        # Add grid point and animation legend
        ax.legend(loc="upper left")
        plt.title(f"Gillespie simulation on grid. {n_rows} x {n_cols}")
        plt.show()
    

     
    #save animation
    writervideo = animation.FFMpegWriter(fps=10) 
    ani.save('gillespie_simulation.mp4', writer=writervideo) 
    plt.close() 

# Visualize the running of the Markov chain on the grid (6x8 grid)
visualize_reaction(trajectory, n_rows=n*2, n_cols=n*2, dynamic_playback=False)



def compute_occupancy_distribution(trajectory, lattice_size, initial_state, time):
    """
    Takes a trajectory of a process on a lattice and computes where the systems spends time
    """

    counter = np.zeros(lattice_size**2)
    time_spent = np.diff(time)
    print(f"{time}")
    idx = 0 #the time_spent array should be traversed sequentially, but the trajectory array should not

    for node in trajectory:
        counter[node] += time_spent[idx]
        print(f"{counter[node]=} and {time_spent[idx]=}")
        idx += 1



    n_steps = len(trajectory)
    plt.bar(np.arange(lattice_size**2), counter, edgecolor='black')
    plt.xticks(np.arange(lattice_size**2), labels=np.arange(lattice_size**2))
    plt.title(f"distribution of state occupancy. {n_steps=}, {gamma_int=}, {gamma_ext=}")
    plt.ylabel("counts")
    plt.xlabel("states on the grid")
    #plt.show()
    
    plt.clf()
    reshaped_counter = counter.reshape((lattice_size,lattice_size))
    plt.imshow(reshaped_counter)
    plt.colorbar()
    plt.title(f"heatmap of occupancy for a {n}x{n} lattice. {n_steps=}, {gamma_int=}, {gamma_ext=}, {initial_state=}")
    plt.savefig(f"gillespie_occupancy_{n}by{n}.png")
    plt.show()
    
    return counter

compute_occupancy_distribution(trajectory, lattice_size=2*n, initial_state=initial_condition, time=time)
