import math
import random
import Speedpos_to_orbit_param as orb_param
import matplotlib.pyplot as plt

constG = 6.67408e-11
UA = 149597900000.0
intervalle = 86400 * 2
annee = 31536000
file_name = 'positions.save'
nb_objets = 600
iterations = 50

coordx_pltl, coordy_pltl = [], []
a = 149597900000.0 * 1
lim = a * 3
with open(file_name, 'r') as file:  # A ENLEVER POUR LE CODE FINAL
    file_lines = file.readlines()
orbit_return_posx, orbit_return_posy = [], []


def orbits(number, i):

    global orbit_return_posx, orbit_return_posy
    global orbit
    orbit = []
    orbit_nu, orbit_a, orbit_e = [], [], []
    orbit_tolerance = 1e-14
    compte = 0
    iterat = i
    while compte < number:

        j = iterat * 4 + 5
        orbit_posx = float(file_lines[j])
        orbit_posy = float(file_lines[j + 1])
        orbit_vx = float(file_lines[j + 2])
        orbit_vy = float(file_lines[j + 3])
        orb_return = orb_param.calculate((orbit_posx, orbit_posy, 0), (orbit_vx, orbit_vy, 0), True)
        orbit_a = abs(float(orb_return[0]))
        orbit_e = abs(float(orb_return[1]))

        # orbital param initialisation
        orbit_E = math.fmod(random.uniform(0, 200 * math.pi), 2 * math.pi)
        orbit_dE = orbit_tolerance + 1
        # orbit_E = orbit_M
        """orbit_E0 = orbit_M
        while orbit_dE > orbit_tolerance:
            orbit_E = orbit_E - (orbit_E - orbit_e * math.sin(orbit_E) - orbit_M) / \
                      (1 - orbit_e * math.cos(orbit_E))
            orbit_dE = math.fabs(orbit_E - orbit_E0)
            orbit_E0 = orbit_E"""
        # if orbit_e >=1: break
        orbit_nu = 2 * math.atan2(math.sqrt(1 + orbit_e) * math.sin(orbit_E / 2),
                                  math.sqrt(abs(1 - orbit_e)) * math.cos(orbit_E / 2))

        orbit_r = orbit_a * (1 - orbit_e * math.cos(orbit_E))

        orbit_return_posx.append(orbit_r * math.cos(orbit_nu))
        orbit_return_posy.append(orbit_r * math.sin(orbit_nu))
        compte += 1

    print(orbit_e, orbit_a)
    return orbit_return_posx, orbit_return_posy


temp_arr_coordx = []
temp_arr_coordy = []

for l in range(nb_objets):
    temp_arr_coordx, temp_arr_coordy = orbits(iterations, l)

coordx_pltl.append(temp_arr_coordx)
coordy_pltl.append(temp_arr_coordy)
nb_objets = nb_objets + iterations * nb_objets


figure = plt.figure()
axes = figure.add_subplot(111)
plt.grid(True)
plt.plot(coordx_pltl, coordy_pltl, linestyle='None', color='g', marker='.')  # planetesimals coordinates
plt.plot([0], marker='o', color='r')  # star's coordinates
plt.tight_layout()
plt.show()
