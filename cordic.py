from fixedpoint import FixedPoint
from math import log2

corners = [
        0.7853981633974483,
        0.4636476090008061,
        0.2449786631268641,
        0.1243549945467614,
        0.0624188099959574,
        0.0312398334302683,
        0.0156237286204768,
        0.0078123410601011,
        0.0039062301319670,
        0.0019531225164788,
        0.0009765621895593,
        0.0004882812111949,
        0.0002441406201494,
        0.0001220703118937,
        6.10351561742e-05,
        3.05175781155e-05,
        1.52587890613e-05,
        7.6293945311e-06,
        3.8146972656e-06,
        1.9073486328e-06,
        9.536743164e-07,
        4.768371582e-07,
        2.384185791e-07,
        1.192092896e-07,
        5.96046448e-08,
        2.98023224e-08,
        1.49011612e-08,
        7.4505806e-09,
        3.7252903e-09,
        1.8626451e-09,
        9.313226e-10,
        4.656613e-10,
    ]

frac_bits = 4
_format = {'signed': True, 'm': 4, 'n': frac_bits}

pi = 3.1415926535897932
pi_fixed = FixedPoint(pi, **_format)
inf = float('inf')


def angles(qf: dict) -> list:
    res = []

    for i in range(len(corners)):
        r = FixedPoint(corners[i], **qf)
        if r == 0:
            break
        res.append(r)
    return res


def k_scale(iters: int, qf: dict) -> FixedPoint:
    k = FixedPoint(1.0, **qf)

    n_bits = qf['n']
    eps = 2 ** -n_bits
    bord = int(0.5 * (log2(1 - 2 * eps + eps ** 2) - log2(2 - eps) + n_bits)) + 1

    for i in range(min(bord, iters)):
        k *= FixedPoint((1 + 2 ** (-2 * i)) ** -0.5, **qf)
    return k


def cordic(angle: FixedPoint, qf: dict):
    x, y = FixedPoint(1, **qf), FixedPoint(0, **qf)
    rotates = angles(qf)
    z = angle

    i = 0
    for i in range(qf['n'] + 1):
        tan = FixedPoint(2 ** -i, **qf)
        if z == 0 or tan == 0:
            break
        sign = 1 if z > 0 else -1

        new_x = x - sign * tan * y
        new_y = y + sign * tan * x
        z -= sign * rotates[i]

        x, y = new_x, new_y
    k = k_scale(i, qf)
    return x * k, y * k


def cordic_asin(y__: FixedPoint, qf: dict):
    if y__ < -1 or y__ > 1:
        raise ValueError("Недопустимое значение")

    rotates = angles(qf)
    x = k_scale(qf['n'] + 1, qf)
    y = FixedPoint(0, **qf)
    z = FixedPoint(0, **qf)

    for i in range(len(rotates)):
        tan = FixedPoint(2 ** -i, **qf)
        if y == y__ or tan == 0:
            break
        sign = 1 if y < y__ else -1

        new_x = x - sign * tan * y
        new_y = y + sign * tan * x
        z += sign * rotates[i]

        x, y = new_x, new_y

    return z


def convert(angle: float, qf: dict) -> tuple:
    global pi
    k = 1

    while angle > pi / 2:
        angle -= pi
        k *= -1
    angle = FixedPoint(angle, **qf)
    return angle, k


def cos(value: float):
    if value < 0:
        return cos(-value)

    value, k = convert(value % (2 * pi), _format)
    return k * cordic(value, _format)[0]


def sin(value: float):
    if value < 0:
        return -sin(-value)

    value, k = convert(value % (2 * pi), _format)
    return k * cordic(value, _format)[1]


def asin(value: float):
    if value < 0:
        return -asin(-value)
    return cordic_asin(FixedPoint(value, **_format), _format)


def acos(value: float):
    return pi_fixed * 0.5 - asin(value)


print(sin(pi / 6))
