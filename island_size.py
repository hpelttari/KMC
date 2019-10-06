import numpy as np


def is_in_range(x, limit):
    return (x >= 0 and x < limit)


def is_not_visited(visited_sites, i, j):
    return not visited_sites[i][j]


def is_correct_height(grid, i, j, height):
    return grid[i][j]==height


def is_valid_to_visit(grid, i, j, height, visited_sites, row, col):
    return (is_in_range(i, row) and
            is_in_range(j, col) and
            is_not_visited(visited_sites, i, j) and
            is_correct_height(grid, i, j, height))


def depth_first_search(grid, i, j, visited_sites, height, row, col, island_size):

    neighbour_row_and_col_indices = [(1, -1),
                                     (0, -1),
                                     (-1, -1),
                                     (-1, 0),
                                     (-1, 1),
                                     (0, 1),
                                     (1, 1),
                                     (1, 0)]

    visited_sites[i][j] = True

    for row_num, col_num in neighbour_row_and_col_indices:
        row_num += i
        col_num += j
        if is_valid_to_visit(grid, row_num, col_num, height, visited_sites, row, col):
            island_size += 1
            island_size = depth_first_search(grid,
                                             row_num,
                                             col_num,
                                             visited_sites,
                                             height,
                                             row,
                                             col,
                                             island_size)

    return island_size


def count_islands(grid, row, col, height):
    visited_sites = np.zeros((row, col))

    count = 0
    island_size = 0
    sizes = []

    for i in range(row):
        for j in range(col):
            if not visited_sites[i][j] and grid[i][j] == height:
                sizes.append(depth_first_search(grid,
                                                i,
                                                j,
                                                visited_sites,
                                                height,
                                                row,
                                                col,
                                                1))
                count += 1

    return count, sizes


def get_island_count_and_average_size(grid):
    row, col = grid.shape
    max_height = np.max(grid)
    island_count = 0
    island_sizes = []

    for height in range(1, int(max_height+1)):
        count, size = count_islands(grid, row, col, height)
        island_count += count
        island_sizes.extend(size)

    return island_count, sum(island_sizes)/island_count


def get_grid_info(grid):
    island_count, avg_size = get_island_count_and_average_size(grid)
    max_height = np.max(grid)
    min_height = np.min(grid)
    height_difference = max_height-min_height
    adatoms = np.sum(grid)
    return island_count, avg_size, max_height, min_height, height_difference, adatoms


def print_grid_info(grid):
    count, size, max_h, min_h, diff, adatoms = get_grid_info(grid)
    print(f'{count} islands, average island size = {size}')
    print(f'number of adatoms = {adatoms}, min height = {min_h}')
    print(f'max height = {max_h}, max height difference {diff}')
