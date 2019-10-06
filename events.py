import numpy as np
from numpy import random
from scipy import constants as const
from functools import partial

def adhere_to_boundary_conditions(grid, x, y):
    boundary = len(grid)
    if x == boundary:
        x = 0
    if y == boundary:
        y = 0
    return x, y

def adsorption(grid, x, y):
    #print("ads!")
    x, y = adhere_to_boundary_conditions(grid, x, y)
    grid[x][y] = grid[x][y] + 1
    return grid


def desorption(grid, x, y):
    x, y = adhere_to_boundary_conditions(grid, x, y)
    grid[x][y] = grid[x][y] - 1
    return grid


def diffusion(grid, x, y, x_dir, y_dir):
    #print("diff!")
    grid[x][y] -=  1
    x2, y2 = adhere_to_boundary_conditions(grid, x+x_dir, y+y_dir)
    grid[x2][y2] += 1
    return grid


def calculate_adsorption_rate(grid, F):
    p = len(grid)**2 * F
    return p


def calculate_diffusion_rate(grid, x, y, E, T, k0):
    n_neighbours = 0
    for i in [-1, 1]:#range(-1, 2):
        ix, iy = adhere_to_boundary_conditions(grid, x+i, y+i)
        if grid[x][iy] >= grid[x][y]:
            n_neighbours += 1
        if grid[ix][y] >= grid[x][y]:
            n_neighbours += 1

    k = k0*np.exp((-E - n_neighbours*E)/(const.k * T))

    return k

def calculate_diffusion_probability(grid, x, y, x_dir, y_dir, p):
    grid_length = len(grid)

    x2, y2 = adhere_to_boundary_conditions(grid, x+x_dir, y+y_dir)
    # check boundary conditions
    """
    if x+x_dir >= grid_length:
        x_dir = -x
    if y+y_dir >= grid_length:
        y_dir = -y
    """
    # diffusion is only possible to the same or lower level
    # as the diffusing atom
    if grid[x2][y2] < grid[x][y]:
        return p

#    print(grid[x][y], grid[x2][y2], "p=0")
    return 0


def calculate_desorption_probability(grid, x, y, p):

    # Desorption is only possible for the adatoms
    if grid[x][y]>0:
        return p

    return 0


def calculate_rates_for_location(grid, x, y, p_adsorption, p_desorption, E, T, k0):
    """
    Calculates the rates of events for a specific location in the grid.

    The necessary rates to calculate are the diffusion rates
    for all possible directions,
    the adsorption rate and the desorption rate.
    """

    # there are 10 rates because the possible events are adsorption,
    # desorption and diffusion,
    # which can happen to 8 different directions
    rates = np.zeros(6)

    rates[0] = p_adsorption
    rates[1] = calculate_desorption_probability(grid, x, y, p_desorption)

    index = 2
    """
    # calculate rates for diffusion
    for i in range(-1, 2):
        for j in range(-1, 2):
            # diffusion to current position is not possible
            if i == 0 and j == 0:
                continue
            rates[index] = calculate_diffusion_probability(grid, x, y, i, j, E, T, k0)
            index+=1
    """
    p_diffusion = calculate_diffusion_rate(grid, x, y, E, T, k0)
    for i in [-1, 1]:
        rates[index] = calculate_diffusion_probability(grid, x, y, i, 0, p_diffusion)
        rates[index+2] =  calculate_diffusion_probability(grid, x, y, 0, i, p_diffusion)
        index += 1

    return rates


def get_direction_from_event_index(event_index):
    if event_index < 2:
        x = 0
        y = 0
    else:
        current_index = 2
        for i in range(-1, 2):
            for j in range(-1, 2):
                if i == 0 and j== 0:
                    continue
                if event_index == current_index:
                    x = i
                    y = j
                    break
                current_index += 1
    return x, y


def calculate_rates_in_local_enviroment(grid, x, y, p_adsorption, p_desorption, rates_list, E, T, k0):
    grid_length = len(grid)
    for i in range(x-1, x+2):
        for j in range(y-1, y+2):
            if i<grid_length and j<grid_length:
                rates_list[i][j] = calculate_rates_for_location(grid, x, y, p_adsorption, p_desorption, E, T, k0)
                continue
            i_1=i
            j_1=j
            if i_1>=grid_length:
                i_1=0
            if j>=grid_length:
                j_1=0
            rates_list[i_1][j_1] = calculate_rates_for_location(grid, x, y, p_adsorption, p_desorption, E, T, k0)

    return rates_list


def update_rates_list(grid, x, y, event_index, rates_list, p_adsorption, p_desorption, p_diffusion):
    x_dir, y_dir = get_direction_from_event_index(event_index)
    x2 = x + x_dir
    y2 = y + y_dir
    x2, y2 = adhere_to_boundary_conditions(grid, x+x_dir, y+y_dir)
    rates_list = calculate_rates_in_local_enviroment(grid, x, y, p_adsorption, p_desorption, p_diffusion, rates_list)
    if x != x2 or y != y2:
        rates_list = calculate_rates_in_local_enviroment(grid, x2, y2, p_adsorption, p_desorption, p_diffusion, rates_list)

    return rates_list

def get_new_rates(grid, rates_list, p_adsorption, p_desorption, E, T, k0):
    for i in range(len(grid)):
        for j in range(len(grid)):
            rates_list[i][j] = calculate_rates_for_location(grid, i, j, p_adsorption, p_desorption, E, T, k0)

    return rates_list

def get_normalized_rates_list(rates_list):
    flattened_rates = rates_list.flatten()
    A = np.sum(flattened_rates)
    return np.cumsum(flattened_rates)/A, A


def choose_event(rates_list):
    r, A = get_normalized_rates_list(rates_list)
    u = random.random()

    for i in range(1, len(r)):
        if r[i-1] < u and u < r[i]:
            return np.unravel_index(i, rates_list.shape), A

    return np.unravel_index(0, rates_list.shape), A


def realize_event(grid, event_indices, events_list):
    x, y, event = event_indices
    return events_list[event](grid, x, y)


def get_delta_t(A):
    u = random.random()
    dt = -(1/A)*np.log(u)
    return dt


def create_events_list():
    event_list = [adsorption, desorption]
    """
    for i in range(-1, 2):
        for j in range(-1, 2):
            if i == 0 and j == 0:
                continue
            event_list.append(partial(diffusion, x_dir=i, y_dir=j))
    """
    event_list.append(partial(diffusion, x_dir=-1, y_dir=0))
    event_list.append(partial(diffusion, x_dir=1, y_dir=0))
    event_list.append(partial(diffusion, x_dir=0, y_dir=-1))
    event_list.append(partial(diffusion, x_dir=0, y_dir=1))
    return event_list
