import numpy as np


def adsorption(grid, x, y):
    grid[x][y] = grid[x][y] + 1
    return grid


def desorption(grid, x, y):
    grid[x][y] = grid[x][y] - 1
    return grid


def diffusion(grid, x, y, x2, y2):
    grid[x][y] = grid[x][y] - 1
    grid[x2][y2] = grid[x2][y2] + 1
    return grid


def calculate_diffusion_probability(grid, x, y, x2, y2, p):

    # diffusion is only possible to the same or lower level
    # as the diffusing atom
    if grid[x2][y2] < grid[x][y]:
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
    for i in range(x-1, x+1):
        for j in range(y-1, y+1):
            # ignore diffusion to current location
            if i==x and j==y:
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

