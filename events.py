import numpy as np
from numpy import random
from scipy import constants as const
from functools import partial

random.seed(0)
def adhere_to_boundary_conditions(grid, x, y):
    boundary = len(grid)
    if x == boundary:
        x = 0
    if y == boundary:
        y = 0
    return x, y

def adsorption(grid, x, y):
    x, y = adhere_to_boundary_conditions(grid, x, y)
    grid[x][y] = grid[x][y] + 1
    return grid


def desorption(grid, x, y):
    x, y = adhere_to_boundary_conditions(grid, x, y)
    grid[x][y] = grid[x][y] - 1
    return grid


def diffusion(grid, x, y, x_dir, y_dir):
    grid[x][y] -=  1
    x2, y2 = adhere_to_boundary_conditions(grid, x+x_dir, y+y_dir)
    grid[x2][y2] += 1
    return grid


def calculate_adsorption_rate(grid, F):
    p = len(grid)**2 * F
    return p


def calculate_diffusion_rate(grid, x, y, E, T, k0):
    n_neighbours = 0
    for i in [-1, 1]:
        ix, iy = adhere_to_boundary_conditions(grid, x+i, y+i)
        if grid[x][iy] >= grid[x][y]:
            n_neighbours += 1
        if grid[ix][y] >= grid[x][y]:
            n_neighbours += 1

    k = k0*np.exp((-E - n_neighbours*E)/(const.k * T))

    return k


def calculate_diffusion_probability(grid, x, y, x_dir, y_dir, p):

    x2, y2 = adhere_to_boundary_conditions(grid, x+x_dir, y+y_dir)

    # diffusion is only possible to the same or lower level
    # as the diffusing atom
    if grid[x2][y2] < grid[x][y]:
        return p

    # in other cases diffusion is not possible
    # so the probability is 0
    return 0


def calculate_desorption_probability(grid, x, y, p):

    # Desorption is only possible for the adatoms
    if grid[x][y]>0:
        return p

    # in other cases desorption is not possible
    # so the probability is 0
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
    p_diffusion = calculate_diffusion_rate(grid, x, y, E, T, k0)

    for i in [-1, 1]:
        rates[index] = calculate_diffusion_probability(grid,
                                                       x,
                                                       y,
                                                       i,
                                                       0,
                                                       p_diffusion)
        rates[index+2] = calculate_diffusion_probability(grid,
                                                         x,
                                                         y,
                                                         0,
                                                         i,
                                                         p_diffusion)
        index += 1

    return rates


def get_direction_from_event_index(event_index):

    indices = [(0, 0), (0, 0), (-1, 0), (1, 0), (0, -1), (0, 1)]
    return indices[event_index]


def calculate_rates_in_local_enviroment(grid, x, y, p_adsorption, p_desorption, rates_list, E, T, k0):
    grid_length = len(grid)

    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            x_i = x+i
            y_j = y+j
            if x_i >= grid_length:
                x_i = 0
            if y_j >= grid_length:
                y_j = 0
            rates_list[x_i][y_j] = calculate_rates_for_location(grid, x_i, y_j, p_adsorption, p_desorption, E, T, k0)

    return rates_list


def update_rates_list(grid, x, y, event_index, rates_list, p_adsorption, p_desorption, E, T, k0):
    grid_length = len(grid)
    rates_list = calculate_rates_in_local_enviroment(grid,
                                                     x,
                                                     y,
                                                     p_adsorption,
                                                     p_desorption,
                                                     rates_list,
                                                     E,
                                                     T,
                                                     k0)
    x_dir, y_dir = get_direction_from_event_index(event_index)
    x += x_dir
    y += y_dir
    if x >= grid_length:
        x = 0
    if y >= grid_length:
        y = 0
    rates_list = calculate_rates_in_local_enviroment(grid,
                                                     x,
                                                     y,
                                                     p_adsorption,
                                                     p_desorption,
                                                     rates_list,
                                                     E,
                                                     T,
                                                     k0)

    return rates_list

def get_new_rates(grid, rates_list, p_adsorption, p_desorption, E, T, k0):
    for i in range(len(grid)):
        for j in range(len(grid)):
            rates_list[i][j] = calculate_rates_for_location(grid,
                                                            i,
                                                            j,
                                                            p_adsorption,
                                                            p_desorption,
                                                            E,
                                                            T,
                                                            k0)

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
    x_and_y_directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    for x_dir, y_dir in x_and_y_directions:
        event_list.append(partial(diffusion, x_dir=x_dir, y_dir=y_dir))

    return event_list
