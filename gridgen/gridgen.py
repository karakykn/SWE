import numpy as np

"""
this code converts blockmesh to gridgen file. works only for 2d
does not handle boundaries.
"""

file_path = 'polyMesh/'  # Replace with the actual file path

x = np.array([])
y = np.array([])
z = np.array([])

with open(file_path + 'points', 'r') as file:
    lines = file.readlines()
    # Start reading after the header and the first integer
    start_reading = False

    for line in lines:
        # Strip leading and trailing whitespace
        stripped_line = line.strip()

        # Check if the line starts with '(' and contains values
        if stripped_line.startswith('('):
            parts = stripped_line.strip('()\n').split()
            if len(parts) == 3:  # Ensure there are three values in the line
                x = np.append(x, float(parts[0]))
                y = np.append(y, float(parts[1]))
                z = np.append(z, float(parts[2]))

precision = 1e-8
N = int(x.size / 2)
x = x[:N]
y = y[:N]
x = np.round(x / precision) * precision
y = np.round(y / precision) * precision
xS = np.unique(x).size
yS = np.unique(y).size
x_grid = x.reshape((xS, yS))
y_grid = y.reshape((xS, yS))
x_grid_transposed = x_grid.T
y_grid_transposed = y_grid.T
xGrd = x_grid_transposed.flatten()
yGrd = y_grid_transposed.flatten()
zGrd = np.zeros(N)

with open(file_path + 'mesh.grd', 'w') as file:
    # Write the single number on the first line
    file.write(f"{1}\n")

    # Write the three numbers on the second line
    file.write(f"{xS} {yS} {1}\n")

    # Write the arrays
    file.write(' '.join(map(str, xGrd)) + '\n')
    file.write(' '.join(map(str, zGrd)) + '\n')
    file.write(' '.join(map(str, yGrd)) + '\n')