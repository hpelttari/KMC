
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
