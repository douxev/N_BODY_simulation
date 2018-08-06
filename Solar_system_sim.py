import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D

nb_annees = int(input("Number of years : "))
intervalle_nb_jours = 2
# constants :
constG = 6.67408e-11
nbiter = 31536000 * nb_annees
UA = 149597900000.0
intervalle = 86400 * intervalle_nb_jours
lim = 35 * UA


# Solar System sim


class BaryCentre:

    def __init__(self):
        self.sun_mass = 1.9884e+30

        self.mass = mercury.mass + venus.mass + earth.mass + mars.mass + jupiter.mass + saturn.mass + \
                    uranus.mass + neptune.mass + self.sun_mass

        self.dt = intervalle
        self.t = 0

    def __iter__(self):
        self.t_abs = 0
        return self

    def __next__(self):
        self.posx = (mercury.mass * mercury.posx + venus.mass * venus.posx + earth.mass * earth.posx + mars.mass
                     * mars.posx + jupiter.mass * jupiter.posx + saturn.mass * saturn.posx + uranus.mass * uranus.posx
                     + neptune.mass * neptune.posx) / self.mass
        self.posy = (mercury.mass * mercury.posy + venus.mass * venus.posy + earth.mass * earth.posy + mars.mass
                     * mars.posy + jupiter.mass * jupiter.posy + saturn.mass * saturn.posy + uranus.mass * uranus.posy
                     + neptune.mass * neptune.posy) / self.mass
        self.posz = (mercury.mass * mercury.posz + venus.mass * venus.posz + earth.mass * earth.posz + mars.mass
                     * mars.posz + jupiter.mass * jupiter.posz + saturn.mass * saturn.posz + uranus.mass * uranus.posz
                     + neptune.mass * neptune.posz) / self.mass
        self.t_abs += self.dt
        self.t += self.dt
        return self


class Planete:
    """represents a planet"""

    star_mass = 1.9884e+30  # Sun mass
    mu = constG * star_mass

    def __init__(self, name, mass, semimajoraxis, e, inclinaison, mean, long_asc, periap):
        self.name = name
        self.tolerance = 1e-20
        self.dt = intervalle
        self.mass = mass
        self.a = semimajoraxis * UA
        self.e = e
        self.i = inclinaison
        self.mean_anomaly = math.pi / 180 * mean
        self.OMEGA = math.pi / 180 * long_asc
        self.omega = math.pi / 180 * periap

    def __iter__(self):
        self.dt = intervalle
        # orbital param initialisation
        self.dE = self.tolerance + 1
        self.M = math.fmod(self.mean_anomaly + self.dt * math.sqrt(self.mu / self.a ** 3), 2 * math.pi)
        self.E = self.E0 = self.M
        while self.dE > self.tolerance:
            self.E = self.E - (self.E - self.e * math.sin(self.E) - self.M) / (1 - self.e * math.cos(self.E))
            self.dE = math.fabs(self.E - self.E0)
            self.E0 = self.E

        self.r = self.a * (1 - self.e * math.cos(self.E))
        self.nu = 2 * math.atan2(math.sqrt(1 + self.e) * math.sin(self.E / 2), math.sqrt(1 - self.e)
                                 * math.cos(self.E / 2))
        self.posx0 = self.r * math.cos(self.nu)
        self.posy0 = self.r * math.sin(self.nu)
        self.posz0 = 0
        self.vx0 = -math.sin(self.E) * math.sqrt(self.mu * self.a) / self.r
        self.vy0 = math.sqrt(1 - self.e ** 2) * math.cos(self.E) * math.sqrt(self.mu * self.a) / self.r
        self.vz0 = 0

        self.posx = self.posx0 * (
                    math.cos(self.omega) * math.cos(self.OMEGA) - math.sin(self.omega) * math.cos(self.i) * math.sin(
                self.OMEGA)) - self.posy0 * (
                                math.sin(self.omega) * math.cos(self.OMEGA) + math.cos(self.omega) * math.cos(
                            self.i) * math.sin(self.OMEGA))

        self.posy = self.posx0 * (
                    math.cos(self.omega) * math.sin(self.OMEGA) + math.sin(self.omega) * math.cos(self.i) * math.cos(
                self.OMEGA)) + self.posy0 * (math.cos(self.omega) * math.cos(self.i) * math.cos(self.OMEGA) - math.sin(
            self.omega) * math.sin(self.OMEGA))

        self.posz = self.posx0 * math.sin(self.omega) * math.sin(self.i) + self.posy0 * math.cos(self.omega) * math.sin(
            self.i)

        self.vx = self.vx0 * (
                    math.cos(self.omega) * math.cos(self.OMEGA) - math.sin(self.omega) * math.cos(self.i) * math.sin(
                self.OMEGA)) - self.vy0 * (
                              math.sin(self.omega) * math.cos(self.OMEGA) + math.cos(self.omega) * math.cos(
                          self.i) * math.sin(self.OMEGA))

        self.vy = self.vx0 * (
                    math.cos(self.omega) * math.sin(self.OMEGA) + math.sin(self.omega) * math.cos(self.i) * math.cos(
                self.OMEGA)) + self.vy0 * (math.cos(self.omega) * math.cos(self.i) * math.cos(self.OMEGA) - math.sin(
            self.omega) * math.sin(self.OMEGA))

        self.vz = self.vx0 * math.sin(self.omega) * math.sin(self.i) + self.vy0 * math.cos(self.omega) * math.sin(
            self.i)

        # end of orbital initialisation param
        return self

    def accel(self, posx, posy, posz):
        # acceleration
        return [((barycentre.posx - posx) * constG * barycentre.mass / (math.sqrt((barycentre.posx - posx) ** 2 +
                                                                                  (barycentre.posy - posy) ** 2 +
                                                                                  (barycentre.posz - posz) ** 2) ** 3)),
                ((barycentre.posy - posy) * constG * barycentre.mass / (math.sqrt((barycentre.posx - posx) ** 2 +
                                                                                  (barycentre.posy - posy) ** 2 +
                                                                                  (barycentre.posz - posz) ** 2) ** 3)),
                ((barycentre.posz - posz) * constG * barycentre.mass / (math.sqrt((barycentre.posx - posx) ** 2 +
                                                                                  (barycentre.posy - posy) ** 2 +
                                                                                  (barycentre.posz - posz) ** 2) ** 3))]

    def __next__(self):
        # k1
        self.k1s = self.accel(self.posx, self.posy, self.posz)

        self.k1px = self.vx
        self.k1py = self.vy
        self.k1pz = self.vz

        # k2
        self.k2s = self.accel(self.posx + self.k1px * self.dt / 2, self.posy + self.k1py * self.dt / 2,
                              self.posz + self.k1pz * self.dt / 2)

        self.k2px = self.vx + self.k1s[0] * self.dt / 2
        self.k2py = self.vy + self.k1s[1] * self.dt / 2
        self.k2pz = self.vz + self.k1s[2] * self.dt / 2

        # k3
        self.k3s = self.accel(self.posx + self.k2px * self.dt / 2, self.posy + self.k2py * self.dt / 2,
                              self.posz + self.k2pz * self.dt / 2)

        self.k3px = self.vx + self.k2s[0] * self.dt / 2
        self.k3py = self.vy + self.k2s[1] * self.dt / 2
        self.k3pz = self.vz + self.k2s[2] * self.dt / 2

        # k4
        self.k4s = self.accel(self.posx + self.k3px * self.dt, self.posy + self.k3py * self.dt,
                              self.posz + self.k3pz * self.dt)

        self.k4px = self.vx + self.k3s[0] * self.dt
        self.k4py = self.vy + self.k3s[1] * self.dt
        self.k4pz = self.vz + self.k3s[2] * self.dt

        # calculating speed at time + dt
        self.vx = self.vx + self.dt * (self.k1s[0] + 2 * self.k2s[0] + 2 * self.k3s[0] + self.k4s[0]) / 6
        self.vy = self.vy + self.dt * (self.k1s[1] + 2 * self.k2s[1] + 2 * self.k3s[1] + self.k4s[1]) / 6
        self.vz = self.vz + self.dt * (self.k1s[2] + 2 * self.k2s[2] + 2 * self.k3s[2] + self.k4s[2]) / 6
        # calculating position at time + dt
        self.posx = self.posx + self.dt * (self.k1px + 2 * self.k2px + 2 * self.k3px + self.k4px) / 6
        self.posy = self.posy + self.dt * (self.k1py + 2 * self.k2py + 2 * self.k3py + self.k4py) / 6
        self.posz = self.posz + self.dt * (self.k1pz + 2 * self.k2pz + 2 * self.k3pz + self.k4pz) / 6
        return self


# Parameters of each Planet ('Name", MASS, A, e, INCLINATION, mean_anomaly, long_ascending node, arg periapsis)
mercury = Planete('Mercury', 330e+21, 0.38710, 0.205631, 7.0049, 174.796, 48.331, 29.124)
venus = Planete('Venus', 4871e+21, 0.72333, 0.006773, 3.3947, 50.115, 76.68069, 54.85229)
earth = Planete('Earth', 5974e+21, 1, 0.016710, 0, 357.51716, 348.73936, 114.20783)
mars = Planete('Mars', 641e+21, 1.52366, 0.093412, 1.8506, 19.3564, 49.562, 286.537)
jupiter = Planete('Jupiter', 1.899e+27, 5.20336, 0.048393, 1.3053, 18.818, 100.492, 275.066)
saturn = Planete('Saturn', 5.68e+26, 9.53707, 0.054151, 2.4845, 320.346750, 113.642811, 336.013862)
uranus = Planete('Uranus', 8.676e+25, 19.1913, 0.047168, 0.7699, 142.955717, 73.989821, 96.541318)
neptune = Planete('Neptune', 1.03e+26, 30.0690, 0.248808, 17.1417, 267.767281, 131.794310, 265.646853)
barycentre = BaryCentre()

mvt_mercury = iter(mercury)
mvt_venus = iter(venus)
mvt_earth = iter(earth)
mvt_mars = iter(mars)
mvt_jupiter = iter(jupiter)
mvt_saturn = iter(saturn)
mvt_uranus = iter(uranus)
mvt_neptune = iter(neptune)
next_barycentre = iter(barycentre)

mercury_coordx, mercury_coordy, mercury_coordz = [], [], []
venus_coordx, venus_coordy, venus_coordz = [], [], []
earth_coordx, earth_coordy, earth_coordz = [], [], []
mars_coordx, mars_coordy, mars_coordz = [], [], []
jupiter_coordx, jupiter_coordy, jupiter_coordz = [], [], []
saturn_coordx, saturn_coordy, saturn_coordz = [], [], []
uranus_coordx, uranus_coordy, uranus_coordz = [], [], []
neptune_coordx, neptune_coordy, neptune_coordz = [], [], []

print("{ans} ans se sont écoulés.".format(ans=int(barycentre.t_abs / 31536000)))

while barycentre.t <= nbiter:  # SIMULATION LOOP

    next(next_barycentre)
    next(mvt_mercury)
    next(mvt_venus)
    next(mvt_earth)
    next(mvt_mars)
    next(mvt_jupiter)
    next(mvt_saturn)
    next(mvt_uranus)
    next(mvt_neptune)

    mercury_coordx.append(mercury.posx)
    mercury_coordy.append(mercury.posy)
    mercury_coordz.append(mercury.posz)
    venus_coordx.append(venus.posx)
    venus_coordy.append(venus.posy)
    venus_coordz.append(venus.posz)
    earth_coordx.append(earth.posx)
    earth_coordy.append(earth.posy)
    earth_coordz.append(earth.posz)
    mars_coordx.append(mars.posx)
    mars_coordy.append(mars.posy)
    mars_coordz.append(mars.posz)
    jupiter_coordx.append(jupiter.posx)
    jupiter_coordy.append(jupiter.posy)
    jupiter_coordz.append(jupiter.posz)
    saturn_coordx.append(saturn.posx)
    saturn_coordy.append(saturn.posy)
    saturn_coordz.append(saturn.posz)
    uranus_coordx.append(uranus.posx)
    uranus_coordy.append(uranus.posy)
    uranus_coordz.append(uranus.posz)
    neptune_coordx.append(neptune.posx)
    neptune_coordy.append(neptune.posy)
    neptune_coordz.append(neptune.posz)

print("{ans} ans se sont écoulés.".format(ans=int(barycentre.t_abs / 31536000)))

figure = plt.figure()
axes = figure.add_subplot(111, projection='3d')
axes.plot(xs=mercury_coordx, ys=mercury_coordy, zs=mercury_coordz, c='r', linestyle=':')
axes.plot(xs=venus_coordx, ys=venus_coordy, zs=venus_coordz, c='b', linestyle=':')
axes.plot(xs=earth_coordx, ys=earth_coordy, zs=earth_coordz, c='g', linestyle=':')
axes.plot(xs=mars_coordx, ys=mars_coordy, zs=mars_coordz, c='m', linestyle=':')
axes.plot(xs=jupiter_coordx, ys=jupiter_coordy, zs=jupiter_coordz, c='c', linestyle=':')
axes.plot(xs=saturn_coordx, ys=saturn_coordy, zs=saturn_coordz, c='y', linestyle=':')
axes.plot(xs=uranus_coordx, ys=uranus_coordy, zs=uranus_coordz, c='k', linestyle=':')
axes.plot(xs=neptune_coordx, ys=neptune_coordy, zs=neptune_coordz, c='r', linestyle=':')
axes.plot([0], [0], [0], marker='o', color='r')  # star's coordinates
axes.set_xlim(-lim, lim)
axes.set_ylim(-lim, lim)
axes.set_zlim(-lim, lim)
plt.grid(True)
plt.show()
