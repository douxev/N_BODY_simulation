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
iterations = 20

coordx_pltl, coordy_pltl = [], []
lim = 2

with open(file_name, 'r') as file:
    file_lines = file.readlines()
orbit_return_posx, orbit_return_posy = [], []
temps = int(int(file_lines[nb_objets * 4 + 5]) / annee)

def orbits(number, i):

    global orbit_return_posx, orbit_return_posy
    global orbit
    orbit = []
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
        orbit_nu = 2 * math.atan2(math.sqrt(1 + orbit_e) * math.sin(orbit_E / 2),
                                  math.sqrt(abs(1 - orbit_e)) * math.cos(orbit_E / 2))

        orbit_r = orbit_a * (1 - orbit_e * math.cos(orbit_E))

        orbit_return_posx.append(orbit_r * math.cos(orbit_nu) / UA)
        orbit_return_posy.append(orbit_r * math.sin(orbit_nu) / UA)
        compte += 1

    return orbit_return_posx, orbit_return_posy


temp_arr_coordx = []
temp_arr_coordy = []

for l in range(nb_objets):
    temp_arr_coordx, temp_arr_coordy = orbits(iterations, l)
coordx_pltl.append(temp_arr_coordx)
coordy_pltl.append(temp_arr_coordy)


figure = plt.figure()
axes = figure.add_subplot(111)
axes.set_xlim(-lim, lim)
axes.set_ylim(-lim, lim)
plt.grid(True)
plt.plot(coordx_pltl, coordy_pltl, linestyle='None', color='g', marker='.')  # planetesimals coordinates
plt.plot([0], marker='o', color='r')  # star's coordinates
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.title("Simulation after {}years\n{}orbits repopulated with {}planetesimals.".
          format(temps, nb_objets, nb_objets * iterations + nb_objets))
plt.xlabel("Relative position to massive object [AU]")
plt.draw()
plt.show()
