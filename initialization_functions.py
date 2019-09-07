import numpy as np
import events


def create_grid(size):
    grid = np.zeros((size, size))
    return grid


def create_events_dictionary():
    events_dict = {"adsorption": events.adsorption, "desorption": events.desorption, "diffusion": events.diffusion}
    return events_dict


def create_rates_list(grid, adsorption_rate):
    """
    Creates the rates list for the initial empty grid.

    In the beginning with the empty grid, only adsorption
    can happen so other rates are 0. The rates list is formed
    so that it is traversed column by column. Each position in
    the grid has 10 possible events. The first rate is adsorption,
    second is desorption and the following 8 are diffusion.
    As an example, the rate for adsorption for position (x,y)
    can be found on index (x*grid_size + y)*10. The grid is
    assumed to always be square so the grid size is the length
    of the side of the grid.
    """
    x_max, y_max = grid.shape

    # every point has 10 possible events
    # which determines the size of the
    # list in combination with
    # the size of the grid
    rates_list = np.zeros(x_max*y_max*10)

    # Desorption and adsorption are not possible
    # in the beginning so only adsorption rates are !=0
    for i in range(0, len(rates_list), 10):
        rates_list[i] = adsorption_rate

    return rates_list
