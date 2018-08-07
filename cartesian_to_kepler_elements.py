import math

mass_central_obj = 1.9884e+30
constG = 6.67408e-11
mu = constG * mass_central_obj

# module that calculates the orbit parameters from speed and pos arrays


def vectorprod(tabl1, tabl2):
    # Vector product
    return (tabl1[1] * tabl2[2] - tabl1[2] * tabl2[1],
            tabl1[2] * tabl2[0] - tabl1[0] * tabl2[2],
            tabl1[0] * tabl2[1] - tabl1[1] * tabl2[0])


def norm(tabl):
    # Norm of vectors
    return math.sqrt(tabl[0] ** 2 + tabl[1] ** 2 + tabl[2] ** 2)


def scalar(tabl1, tabl2):
    # Scalar product
    return tabl1[0] * tabl2[0] + tabl1[1] * tabl2[1] + tabl1[2] * tabl2[2]


def calculate(position, speed, two_dim):
    # to calculate parameters
    h = vectorprod(position, speed)
    e = vectorprod(speed, h)[0] / mu - position[0] / norm(position), \
        vectorprod(speed, h)[1] / mu - position[1] / norm(position), \
        vectorprod(speed, h)[2] / mu - position[2] / norm(position)

    n = (-h[1], h[0], 0)

    if scalar(position, speed) >= 0:
        nu = math.acos(scalar(e, position) / norm(e) / norm(position))
    else:
        nu = 2 * math.pi - math.acos(scalar(e, position) / norm(e) / norm(position))

    var_e = norm(e)
    E = 2 * math.atan(math.tan(nu / 2) / math.sqrt((1+var_e) / (1 - var_e)))
    mean = E - var_e * math.sin(E)
    a = 1 / (2 / norm(position) - norm(position) ** 2 / mu)

    if not two_dim:
        i = math.acos(h[2] / norm(h))
        if n[1] >= 0:
            OMEGA = math.acos(n[0] / norm(n))
        else:
            OMEGA = 2 * math.pi - math.acos(n[0] / norm(n))

        if e[2] >= 0:
            omega = math.acos(scalar(n, e) / norm(n) / norm(e))
        else:
            omega = 2 * math.pi - math.acos(scalar(n, e) / norm(n) / norm(e))


        return a, var_e, mean, i, omega, OMEGA
    else:
        return float(a), float(var_e), float(mean)

# print("\nSemi-major axis a = {a}\nEccentricity e = {e}\nArgument of periapsis omega = {omega}\nLAN OMEGA = {OMEGA}\n"
#       "Inclination i = {i}\nMean anomaly M = {mean}\n".format(a=a, e=var_e, omega=omega, OMEGA=OMEGA, i=i, mean=mean))
