from fixedpoint import FixedPoint
from math import log2, sqrt
from matplotlib import pyplot
from time import time
import numpy as np

from math import sin as math_sin
from math import cos as math_cos
from math import asin as math_asin

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

frac_bits = 12
_format = {'signed': True, 'm': 8, 'n': frac_bits, 'rounding': 'in'}

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
    tan = FixedPoint(1.0, **qf)
    for i in range(len(rotates)):
        if z == 0 or tan == 0:
            break
        sign = 1 if z > 0 else -1

        new_x = x - sign * tan * y
        new_y = y + sign * tan * x
        z -= sign * rotates[i]

        x, y = new_x, new_y
        tan >>= 1

    k = k_scale(i, qf)
    return k * x, k * y


def cordic_asin(y__: FixedPoint, qf: dict):
    if y__ < -1 or y__ > 1:
        raise ValueError("Недопустимое значение")

    rotates = angles(qf)
    x = k_scale(qf['n'] + 1, qf)
    y = FixedPoint(0, **qf)
    z = FixedPoint(0, **qf)

    i = 0
    tan = FixedPoint(1.0, **qf)
    for i in range(len(rotates)):
        if y == y__ or tan == 0:
            break
        sign = 1 if y < y__ else -1

        new_x = x - sign * tan * y
        new_y = y + sign * tan * x
        z += sign * rotates[i]

        x, y = new_x, new_y
        tan >>= 1

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
    return float(k * cordic(value, _format)[0])


def sin(value: float):
    if value < 0:
        return -sin(-value)

    value, k = convert(value % (2 * pi), _format)
    return float(k * cordic(value, _format)[1])


def asin(value: float):
    if value < 0:
        return -asin(-value)
    return float(cordic_asin(FixedPoint(value, **_format), _format))


def acos(value: float):
    return float(pi_fixed * 0.5 - asin(value))


def iterations():
    pyplot.title("Время работы в зависимости от количества бит ")
    pyplot.grid()
    pyplot.xlabel('Количество бит')
    pyplot.ylabel('Количество итераций')

    format_n = _format['n']

    x = tuple(range(1, 25))
    y = tuple()

    for i in x:
        _format['n'] = i

        s = time()
        res = cordic_asin(FixedPoint(0.423), _format)
        y += (time() - s,)

    # pyplot.plot((x[0], x[-1]), (y[0], y[-1]), color='red')
    pyplot.plot((x[0], x[-1]), (y[0], y[-1]), color='red')
    pyplot.plot(x, y, color='black')
    pyplot.show()

    _format['n'] = format_n


def graphics(
        x_0: float,
        x_n: float,
        cnt: int,
        cordic_func,
        math_func,
        title_form: str,
        file_form: str,
        x_label: str,
        func_title: str,
        list_fracs: [list | tuple] = tuple(range(1, 9))
):
    format_n = _format['n']

    x = np.linspace(x_0, x_n, cnt)
    y_math = tuple(map(lambda i: math_func(i), x))

    pyplot.xlabel(x_label)
    pyplot.ylabel('value')

    for n in list_fracs:
        _format['n'] = n
        y_cordic = tuple(map(lambda i: cordic_func(i), x))

        pyplot.title(title_form.format(n))
        pyplot.xlabel(x_label)
        pyplot.ylabel('value')
        pyplot.grid()

        pyplot.plot(x, y_math, linestyle='--', color='black', label=func_title)
        pyplot.plot(x, y_cordic, color='red', label=f'{func_title}_cordic')
        pyplot.legend()

        pyplot.savefig(file_form.format(n))

        pyplot.legend().remove()
        pyplot.cla()

    pyplot.clf()
    _format['n'] = format_n


tests = (dict(x_0=-pi / 2,
              x_n=pi / 2,
              cnt=100,
              cordic_func=sin,
              math_func=math_sin,
              title_form='sin {} bits',
              file_form='Graphics\\sin\\sin_{}.png',
              x_label='radians',
              func_title='sin'),
         dict(x_0=0,
              x_n=pi,
              cnt=100,
              cordic_func=cos,
              math_func=math_cos,
              title_form='cos {} bits',
              file_form='Graphics\\cos\\cos_{}.png',
              x_label='radians',
              func_title='cos'),
         dict(x_0=-1.0,
              x_n=1.0,
              cnt=100,
              cordic_func=asin,
              math_func=math_asin,
              title_form='arcsine {} bits',
              file_form='Graphics\\arcsine\\arcsine_{}.png',
              x_label='number',
              func_title='arcsine')
         )


def main():
    for info in tests:
        graphics(**info)


main()