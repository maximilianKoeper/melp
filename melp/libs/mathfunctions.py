import numpy as np


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2' """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return (np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) / np.pi) * 180


def distance_between_2d(p1, p2):
    """ Return the distance between two points """
    return np.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)


def distance_between_3d(p1, p2):
    """ Return the distance between two points """
    return np.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2)


def angle_between_phi(v1, v2):
    """ Returns the phi angle between two vectors """
    v1 = unit_vector(v1)
    v2 = unit_vector(v2)
    tmp = np.arctan2(v1[0] * v2[1] - v1[1] * v2[0], v1[0] * v2[0] + v1[1] * v2[1])

    return (tmp / np.pi) * 180
