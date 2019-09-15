import numpy as np
from numpy import random
from functools import partial

def adsorption(grid, x, y):
    grid[x][y] = grid[x][y] + 1
    return grid


def desorption(grid, x, y):
    grid[x][y] = grid[x][y] - 1
    return grid


def diffusion(grid, x, y, x_dir, y_dir):
    grid[x][y] -=  1
    grid[x+x_dir][y+y_dir] += 1
    return grid


def calculate_diffusion_probability(grid, x, y, x_dir, y_dir, p):
    grid_length = len(grid)

    # check boundary conditions
    if x+x_dir >= grid_length:
        x_dir = -x
    if y+y_dir >= grid_length:
        y_dir = -y

    # diffusion is only possible to the same or lower level
    # as the diffusing atom
    if grid[x+x_dir][y+y_dir] < grid[x][y]:
        return p

    return 0


def calculate_desorption_probability(grid, x, y, p):

    # Desorption is only possible for the adatoms
    if grid[x][y]>0:
        return p

    return 0


def calculate_rates_for_location(grid, x, y, p_adsorption, p_desorption, p_diffusion):
    """
    Calculates the rates of events for a specific location in the grid.

    The necessary rates to calculate are the diffusion rates
    for all possible directions,
    the adsorption rate and the desorption rate.
    """

    # there are 10 rates because the possible events are adsorption,
    # desorption and diffusion,
    # which can happen to 8 different directions
    rates = np.zeros(10)

    rates[0] = p_adsorption
    rates[1] = calculate_desorption_probability(grid, x, y, p_desorption)

    index = 2
    # calculate rates for diffusion
    for i in range(-1, 2):
        for j in range(-1, 2):
            # diffusion to current position is not possible
            if i == j:
                continue
            rates[index] = calculate_diffusion_probability(grid, x, y, i, j, p_diffusion)
            index+=1

    return rates


def calculate_rates_in_local_enviroment(grid, x, y, p_adsorption, p_desorption, p_diffusion, rates_list):
    grid_length = len(grid)
    for i in range(x-1, x+2):
        for j in range(y-1, y+2):
            if i<grid_length and j<grid_length:
                rates_list[i][j] = calculate_rates_for_location(grid, x, y, p_adsorption, p_desorption, p_diffusion)
                continue
            i_1=i
            j_1=j
            if i_1>=grid_length:
                i_1=0
            if j>=grid_length:
                j_1=0
            rates_list[i_1][j_1] = calculate_rates_for_location(grid, x, y, p_adsorption, p_desorption, p_diffusion)

    return rates_list


def update_rates_list(grid, x, y, x2, y2, rates_list, p_adsorption, p_desorption, p_diffusion):

    rates_list = calculate_rates_in_local_enviroment(grid, x, y, p_adsorption, p_desorption, p_diffusion, rates_list)
    rates_list = calculate_rates_in_local_enviroment(grid, x2, y2, p_adsorption, p_desorption, p_diffusion, rates_list)

    return rates_list


def get_normalized_rates_list(rates_list):
    flattened_rates = rates_list.flatten()
    A = np.sum(flattened_rates)
    return np.cumsum(flattened_rates)/A


def choose_event(rates_list):
    r = get_normalized_rates_list(rates_list)
    u = random.random()

    for i in range(1, len(r)):
        if r[i-1] < u and u < r[i]:
            return np.unravel_index(i, rates_list.shape)

    return np.unravel_index(0, rates_list.shape)


def realize_event(grid, event_indices, events_list):
    x, y, event = event_indices
    return events_list[event](grid, x, y)


def create_events_list():
    event_list = [adsorption, desorption]

    for i in range(-1, 2, 2):
        for j in range(-1, 2, 2):
            event_list.append(partial(diffusion, x_dir=i, y_dir=j))

    return event_list
