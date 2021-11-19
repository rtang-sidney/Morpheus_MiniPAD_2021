import numpy as np

AXIS_X = "x"
AXIS_Y = "y"
AXIS_Z = "z"
VALID_AXES = [AXIS_X, AXIS_Y, AXIS_Z]

d = 3.355e-10
wavelength = 4.905e-10
theta = np.arcsin(wavelength / (2 * d))


def _angles_incoming(axis, delta_i):
    """
    return the two turning angles of the polarisation of the incoming beam
    :param axis: x, y or z, the incoming polarisation in the scattering frame
    :param delta_i: the angle between ki and Q vectors, equal to theta at an elastic scattering
    :return: alpha_i and beta_i, the turning angles around yi and zi axes, respectively
    """
    if axis == AXIS_X:
        alpha_i = np.pi / 2.0
        beta_i = -delta_i
    elif axis == AXIS_Y:
        alpha_i = np.pi / 2.0
        beta_i = np.pi / 2.0 - delta_i
    elif axis == AXIS_Z:
        alpha_i = 0
        beta_i = 0
    else:
        raise ValueError("Invalid axis of the incoming polarisation given: {}".format(axis))
    return alpha_i, beta_i


def _angles_outgoing(axis, delta_f):
    """
    return the two turning angles of the polarisation of the outgoing beam
    :param axis: x, y or z, the outgoing polarisation in the scattering frame
    :param delta_f: the angle between kf and Q vectors, equal to theta at an elastic scattering
    :return: alpha_f and beta_f, the turning angles around yf and zf axes, respectively
    """
    if axis == AXIS_X:
        alpha_f = -np.pi / 2.0
        beta_f = delta_f
    elif axis == AXIS_Y:
        alpha_f = np.pi / 2.0
        beta_f = np.pi / 2.0 + delta_f
    elif axis == AXIS_Z:
        alpha_f = 0
        beta_f = 0
    else:
        raise ValueError("Invalid axis of the outgoing polarisation given: {}".format(axis))
    return alpha_f, beta_f


def _adjust_angle(angle):
    while angle < 0:
        angle += 2 * np.pi
    while angle > 2 * np.pi:
        angle -= 2 * np.pi
    return angle


def polarisation2angles(incoming_axis, outgoing_axis, twotheta):
    """
    return the four turning angles of the polarisation of the incoming and outgoing beams
    :param incoming_axis: x, y or z, the incoming polarisation in the scattering frame
    :param outgoing_axis: x, y or z, the outgoing polarisation in the scattering frame
    :param twotheta: scattering angle at an elastic scattering
    :return: alpha_i, beta_i, alpha_f and beta_f, the turning angles around yi, zi, yf and zf axes, respectively
    """
    if incoming_axis not in VALID_AXES:
        raise ValueError(
            "Invalid axis of the incoming polarisation given: {}. It has to be one of the followings: {}".format(
                incoming_axis, VALID_AXES))
    if outgoing_axis not in VALID_AXES:
        raise ValueError(
            "Invalid axis of the outgoing polarisation given: {}. It has to be one of the followings: {}".format(
                outgoing_axis, VALID_AXES))
    alpha_i, beta_i = _angles_incoming(axis=incoming_axis, delta_i=twotheta / 2.0)
    alpha_f, beta_f = _angles_outgoing(axis=outgoing_axis, delta_f=np.pi - twotheta / 2.0)
    alpha_f = -alpha_f
    beta_f = -beta_f
    return _adjust_angle(alpha_i), _adjust_angle(beta_i), _adjust_angle(alpha_f), _adjust_angle(beta_f)
