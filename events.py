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
    """
    Calculates the rates for the location (x,y) and for the locations aroun it.

    These points are the only points which rates might change after an event, so calculating
    the rates for these points only is sufficient in every Monte Carlo step.
    """
    grid_length = len(grid)

    for i in range(x-1, x+1):
        for j in range(y-1, y+1):
            start_index = i*grid_length + j*10
            rates_to_modify = rates_list[start_index:start_index+10]
            rates = calculate_rates_for_location(grid, x, y, p_adsorption, p_desorption, p_diffusion)
            for index, rate in enumerate(rates):
                rates_to_modify[index] = rate

    return rates_list


def update_rates_list(grid, x, y, x2, y2, rates_list, p_adsorption, p_desorption, p_diffusion):

    rates_list = calculate_rates_in_local_enviroment(grid, x, y, p_adsorption, p_desorption, p_diffusion, rates_list)
    rates_list = calculate_rates_in_local_enviroment(grid, x2, y2, p_adsorption, p_desorption, p_diffusion, rates_list)

    return rates_list








