import matplotlib.pyplot as plt
import math
import random
import os

nb_objets = 600
nb_annees = int(input("Number of years : "))
nb_UA = 1
intervalle_nb_jours = 2
loin = 1.5  # éloignement des planetesimales / planete
proche = 0.5
# excentricity
e = 0
e_pltl = 0
# constants :
constG = 6.67408e-11
nbiter = 31536000 * nb_annees
a = 149597900000.0 * nb_UA
intervalle = 86400 * intervalle_nb_jours
lim = a * 3
file_name = 'positions.save'
reset = int(input("1 to reset the simulation, 0 to continue : "))


# Working on NBODY_RK


class Planetesimal:
    """represents a planetimal"""

    plan_mass = 8.682e+25
    star_mass = 1.9884e+30
    mu = constG * star_mass

    def __init__(self, number):
        self.number = int(number)
        self.tolerance = 1e-14
        self.dt = intervalle

    def __iter__(self):
        self.dt = intervalle
        # orbital param initialisation
        self.a_pltl = random.randrange(a * proche, loin * a, 500)

        self.dE = self.tolerance + 1
        self.M = math.fmod(random.uniform(0, 2 * math.pi) +
                           self.dt * math.sqrt(self.mu / self.a_pltl ** 3), 2 * math.pi)

        self.E = self.E0 = self.M
        while self.dE > self.tolerance:
            self.E = self.E - (self.E - e_pltl * math.sin(self.E) - self.M) / (1 - e_pltl * math.cos(self.E))
            self.dE = math.fabs(self.E - self.E0)
            self.E0 = self.E

        self.nu = 2 * math.atan2(math.sqrt(1 + e_pltl) * math.sin(self.E / 2),
                                 math.sqrt(1 - e_pltl) * math.cos(self.E / 2))

        self.r = self.a_pltl * (1 - e_pltl * math.cos(self.E))

        self.posx = self.r * math.cos(self.nu)
        self.posy = self.r * math.sin(self.nu)

        self.vx = -math.sin(self.E) * math.sqrt(self.mu * self.a_pltl) / self.r
        self.vy = math.sqrt(1 - e_pltl ** 2) * math.cos(self.E) * math.sqrt(self.mu * self.a_pltl) / self.r
        # end of orbital initialisation param
        return self

    def accel_planetesimal(self, posx, posy):
        # acceleration of one planetesimal
        return ((-posx * constG * self.star_mass / (math.sqrt(posx ** 2 + posy ** 2) ** 3)) -
                ((posx - planete.posx) * constG * self.plan_mass / (math.sqrt((planete.posx - posx) ** 2 +
                                                                              (planete.posy - posy) ** 2) ** 3)),
                (-posy * constG * self.star_mass / (math.sqrt(posx ** 2 + posy ** 2) ** 3)) -
                ((posy - planete.posy) * constG * self.plan_mass / (math.sqrt((planete.posx - posx) ** 2 +
                                                                              (planete.posy - posy) ** 2) ** 3)))

    def __next__(self):

        # if (planete.posx - self.posx) ** 2 + (planete.posy - self.posy) ** 2 > 2.7e+15:  # ~50e+6m is global radius
        # k1
        self.k1s = self.accel_planetesimal(self.posx, self.posy)

        self.k1px = self.vx
        self.k1py = self.vy

        # k2
        self.k2s = self.accel_planetesimal(self.posx + self.k1px * self.dt / 2, self.posy + self.k1py * self.dt / 2)

        self.k2px = self.vx + self.k1s[0] * self.dt / 2
        self.k2py = self.vy + self.k1s[1] * self.dt / 2

        # k3
        self.k3s = self.accel_planetesimal(self.posx + self.k2px * self.dt / 2, self.posy + self.k2py * self.dt / 2)

        self.k3px = self.vx + self.k2s[0] * self.dt / 2
        self.k3py = self.vy + self.k2s[1] * self.dt / 2

        # k4
        self.k4s = self.accel_planetesimal(self.posx + self.k3px * self.dt, self.posy + self.k3py * self.dt)

        self.k4px = self.vx + self.k3s[0] * self.dt
        self.k4py = self.vy + self.k3s[1] * self.dt

        # calculating speed at time + dt
        self.vx = self.vx + self.dt * (self.k1s[0] + 2 * self.k2s[0] + 2 * self.k3s[0] + self.k4s[0]) / 6
        self.vy = self.vy + self.dt * (self.k1s[1] + 2 * self.k2s[1] + 2 * self.k3s[1] + self.k4s[1]) / 6
        # calculating position at time + dt
        self.posx = self.posx + self.dt * (self.k1px + 2 * self.k2px + 2 * self.k3px + self.k4px) / 6
        self.posy = self.posy + self.dt * (self.k1py + 2 * self.k2py + 2 * self.k3py + self.k4py) / 6

        return self


class Planete:
    """represents a planet"""

    mass = 8.682e+25  # Uranus' mass
    star_mass = 1.9884e+30  # Sun mass
    mu = constG * star_mass

    def __init__(self, name):
        self.name = name
        self.t = 0
        self.planetesimal_mass = 9.445e+20  # Ceres mass
        self.tolerance = 1e-14
        self.dt = intervalle
        self.t_abs = 0
        self.objects_to_remove = []

    def __iter__(self):
        self.dt = intervalle
        # orbital param initialisation
        self.dE = self.tolerance + 1
        self.M = math.fmod(self.dt * math.sqrt(self.mu / a ** 3), 2 * math.pi)
        self.E = self.E0 = self.M
        while self.dE > self.tolerance:
            self.E = self.E - (self.E - e * math.sin(self.E) - self.M) / (1 - e * math.cos(self.E))
            self.dE = math.fabs(self.E - self.E0)
            self.E0 = self.E

        self.r = a * (1 - e * math.cos(self.E))
        self.nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(self.E / 2), math.sqrt(1 - e) * math.cos(self.E / 2))
        self.posx = a
        self.posy = 0
        self.vx = -math.sin(self.E) * math.sqrt(self.mu * a) / self.r
        self.vy = math.sqrt(1 - e ** 2) * math.cos(self.E) * math.sqrt(self.mu * a) / self.r
        # end of orbital initialisation param
        return self

    def accel(self, posx, posy):
        # acceleration

        self.coord_x_accel = (-posx * constG * self.star_mass / (math.sqrt(posx ** 2 + posy ** 2) ** 3))
        self.coord_y_accel = (-posy * constG * self.star_mass / (math.sqrt(posx ** 2 + posy ** 2) ** 3))

        for k in range(nb_objets):
            self.coord_x_accel += ((planetesimal[k].posx - posx) * constG * self.planetesimal_mass
                                   / (math.sqrt((planetesimal[k].posx - posx) ** 2 +
                                                (planetesimal[k].posy - posy) ** 2) ** 3))
            self.coord_y_accel += ((planetesimal[k].posy - posy) * constG * self.planetesimal_mass
                                   / (math.sqrt((planetesimal[k].posx - posx) ** 2 +
                                                (planetesimal[k].posy - posy) ** 2) ** 3))

        return self.coord_x_accel, self.coord_y_accel

    def __next__(self):
        # k1
        self.k1s = self.accel(self.posx, self.posy)

        self.k1px = self.vx
        self.k1py = self.vy

        # k2
        self.k2s = self.accel(self.posx + self.k1px * self.dt / 2, self.posy + self.k1py * self.dt / 2)

        self.k2px = self.vx + self.k1s[0] * self.dt / 2
        self.k2py = self.vy + self.k1s[1] * self.dt / 2

        # k3
        self.k3s = self.accel(self.posx + self.k2px * self.dt / 2, self.posy + self.k2py * self.dt / 2)

        self.k3px = self.vx + self.k2s[0] * self.dt / 2
        self.k3py = self.vy + self.k2s[1] * self.dt / 2

        # k4
        self.k4s = self.accel(self.posx + self.k3px * self.dt, self.posy + self.k3py * self.dt)

        self.k4px = self.vx + self.k3s[0] * self.dt
        self.k4py = self.vy + self.k3s[1] * self.dt

        # calculating speed at time + dt
        self.vx = self.vx + self.dt * (self.k1s[0] + 2 * self.k2s[0] + 2 * self.k3s[0] + self.k4s[0]) / 6
        self.vy = self.vy + self.dt * (self.k1s[1] + 2 * self.k2s[1] + 2 * self.k3s[1] + self.k4s[1]) / 6
        # calculating position at time + dt
        self.posx = self.posx + self.dt * (self.k1px + 2 * self.k2px + 2 * self.k3px + self.k4px) / 6
        self.posy = self.posy + self.dt * (self.k1py + 2 * self.k2py + 2 * self.k3py + self.k4py) / 6
        self.t_abs += self.dt
        self.t += self.dt
        return self


planete = Planete('Uranus')
print(planete.name)
coordx, coordy = [], []  # lists containing coordinates of the planet

# planetesimal setting
coordx_pltl, coordy_pltl = [], []  # lists containing planetesimal's coordinates
planetesimal = []
mvt_pltl = []

for i in range(nb_objets):
    planetesimal.append(Planetesimal(i))

# FILE READING
with open(file_name, 'r') as file:
    file_lines = file.readlines()
    if os.path.getsize(file_name) > 0 and reset == 0:  # then write the file if file empty
        written = 1
    else:  # if not empty then read previous positions and write new ones
        written = 0

# FILE READING
if written == 1:
    nb_objets = int(file_lines[0])
    planete.posx = float(file_lines[1])
    planete.posy = float(file_lines[2])
    planete.vx = float(file_lines[3])
    planete.vy = float(file_lines[4])

    for i in range(nb_objets):
        j = i * 4 + 5
        planetesimal[i].posx = float(file_lines[j])
        planetesimal[i].posy = float(file_lines[j + 1])
        planetesimal[i].vx = float(file_lines[j + 2])
        planetesimal[i].vy = float(file_lines[j + 3])

        # print("#{i}: posx: {posx}, posy: {posy}\n".format(posx=planetesimal[i].posx,posy=planetesimal[i].posy,i=i))
    planete.t_abs = int(file_lines[nb_objets * 4 + 5])
    intervalle = int(file_lines[nb_objets * 4 + 6])
# ENDING FILE READING
print("{ans} ans se sont écoulés avant le début de la simulation.\n{obj} planétésimales.\n".format(
    ans=int(planete.t_abs / 31536000), obj=nb_objets))

for i in range(nb_objets):  # iter launch
    mvt_pltl.append(iter(planetesimal[i]))
mvt = iter(planete)


while planete.t < nbiter:  # SIMULATION LOOP


    for i in range(nb_objets):
        next(mvt_pltl[i])

    next(mvt)
    # if planete.t > 31536000 * 999000:
    coordx.append(planete.posx)
    coordy.append(planete.posy)


print("{ans} ans se sont écoulés depuis le début du calcul.\n{obj} planétésimales.\n".format(
    ans=int(planete.t_abs / 31536000), obj=nb_objets))
for i in range(nb_objets):
    coordx_pltl.append(planetesimal[i].posx)
    coordy_pltl.append(planetesimal[i].posy)
figure = plt.figure()
axes = figure.add_subplot(111)
axes.set_xlim(-lim, lim)
axes.set_ylim(-lim, lim)
plt.grid(True)
plt.plot(coordx, coordy, linestyle='-.', marker=',')  # planet's coordinates
plt.plot(coordx_pltl, coordy_pltl, linestyle='None', color='g', marker='.')  # planetesimals coordinates
plt.plot([0], marker='o', color='r')  # star's coordinates
plt.tight_layout()
plt.show()

with open(file_name, 'w') as file:
    file.write("{}\n".format(nb_objets))
    file.write("{posx}\n{posy}\n{vx}\n{vy}\n".format(posx=planete.posx, posy=planete.posy,
                                                     vx=planete.vx, vy=planete.vy))
    for i in range(0, nb_objets):
        file.write("{posx2}\n{posy2}\n{vx2}\n{vy2}\n".format(posx2=planetesimal[i].posx, posy2=planetesimal[i].posy,
                                                             vx2=planetesimal[i].vx, vy2=planetesimal[i].vy, ))
    file.write("{t}\n{dt}\n".format(t=planete.t_abs, dt=intervalle))
