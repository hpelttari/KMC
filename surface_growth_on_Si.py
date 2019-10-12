import time

import scipy.constants as const
from matplotlib import pyplot as plt

from initialization_functions import read_config, create_grid, create_rates_list
from events import calculate_adsorption_rate, create_events_list, choose_event, realize_event, get_delta_t, update_rates_list
from island_size import print_grid_info, get_grid_info

config = read_config('config.yml')

island_counts = []
island_sizes = []
max_heights = []
min_heights = []
height_differences = []
num_adatoms = []

for T in config['temperatures']:
    grid = create_grid(config['grid_size'])
    adsorption_rate = calculate_adsorption_rate(grid, config['F'])
    rates = create_rates_list(grid, adsorption_rate)
    events = create_events_list()

    t = 0
    k0 = 2*const.k*T/const.h
    eV = const.physical_constants["electron volt"][0]
    E = config["E"]*eV

    t1 = time.clock()
    while t < 10:
        event_indices, A = choose_event(rates)
        grid = realize_event(grid, event_indices, events)
        t += get_delta_t(A)
        rates = update_rates_list(grid,
                                  *event_indices,
                                  rates,
                                  adsorption_rate,
                                  config['desorption_rate'],
                                  E,
                                  T,
                                  k0)
    print(grid)
    print(f"t = {time.clock()-t1}")

    grid_info = get_grid_info(grid)
    island_counts.append(grid_info[0])
    island_sizes.append(grid_info[1])
    max_heights.append(grid_info[2])
    min_heights.append(grid_info[3])
    height_differences.append(grid_info[4])
    num_adatoms.append(grid_info[5])

    print_grid_info(grid)
    plt.figure(T)
    plt.title(f"Surface growth of Si in {T}K")
    plt.imshow(grid, cmap='gray')
    plt.savefig(f"surface_growth {T}K.jpg")

plt.figure(1)
plt.title("island counts as a function of temperature")
plt.plot(config["temperatures"], island_counts)
plt.savefig("counts.jpg")

plt.figure(2)
plt.title("sizes")
plt.plot(config["temperatures"], island_sizes)
plt.savefig("sizes.jpg")

plt.figure(3)
plt.title("max heigts")
plt.plot(config["temperatures"], max_heights)
plt.savefig("max_heights.jpg")

plt.figure(4)
plt.title("min heights")
plt.plot(config["temperatures"], min_heights)
plt.savefig("min_heights.jpg")

plt.figure(5)
plt.title("height differences")
plt.plot(config["temperatures"], height_differences)
plt.savefig("differences.jpg")

plt.figure(6)
plt.title("number of adatoms")
plt.plot(config["temperatures"], num_adatoms)
plt.savefig("adatoms.jpg")

plt.show()
