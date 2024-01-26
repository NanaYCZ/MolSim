import numpy as np


# smoothed lennard-jones potential implementation in python
# uses no optimizations compared to cpp implementation
# and whenever changes are made to the cpp implementation
# it can be verified with this implementation
def smoothed_LJ(xi, xj, sigma = 1.0, epsilon =1.0, cutoff = 3.0, r_l = 1.5):
    xi = np.array(xi)
    xj = np.array(xj)

    dij = np.linalg.norm(xi - xj)

    lj_force = 0

    if dij <= r_l:
        lj_force = (-24 * epsilon  / dij**2) * ((sigma / dij) ** 6 - 2 * (sigma / dij) ** 12) * (xi - xj)
    elif r_l < dij <= cutoff:
        r_c = cutoff
        prefac = - 24 * (sigma**6 * epsilon * (r_c - dij) / (dij**14 * (r_c - r_l)**3) )
        prefac *= r_c**2 * (2 * sigma**6 - dij**6) + r_c * (3 * r_l - dij) * (dij**6 - 2 * sigma**6) + \
                         dij * (5 * r_l * sigma**6 - 2 * r_l * dij**6 - 3 * sigma**6 * dij + dij**7)

        lj_force = prefac * (xj - xi)
    else:
        # Distance exceeds cutoff, no force
        lj_force = np.zeros(3)

    return lj_force


sigma = 1.0
epsilon = 1.0
cutoff = 3.0
r_l = 1.5
result1 = smoothed_LJ([1.0,1,1],[2.0,1.1,1.2])
result2 = smoothed_LJ([5.7,4.5,1],[5.9,4.0,0.5])
result3 = smoothed_LJ([0.0,0.1,0.01],[0.1,0.7,0.3])

print("1: " + str(result1))
print("2: " + str(result2))
print("3: " + str(result3))


