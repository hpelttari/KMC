
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
