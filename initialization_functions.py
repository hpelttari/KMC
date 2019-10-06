import numpy as np
import yaml

import events


def create_grid(size):
    grid = np.zeros((size, size))
    return grid


def create_events_dictionary():
    events_dict = {"adsorption": events.adsorption, "desorption": events.desorption, "diffusion": events.diffusion}
    return events_dict


def create_rates_list(grid, adsorption_rate):
    grid_size = len(grid)

    rates = np.zeros((grid_size, grid_size, 6))
    for i in range(grid_size):
        for j in range(grid_size):
            rates[i][j][0] = adsorption_rate

    return rates


def read_config(filename):

    with open(filename, 'r') as configuration:
        try:
            config = yaml.safe_load(configuration)
        except yaml.YAMLErro as error:
            print(error)

    return config
