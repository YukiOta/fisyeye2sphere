# -*- coding: utf-8 -*-

import numpy as np


def sphere2cartesian(coord, r=1):
    """球座標をデカルト座標系に変換
    Parameters
    ----------
    coord[0] : 緯度
    coord[1] : 経度

    Retrun
    ----------
    x, y, z : 各座標
    """
    col = coord[0]
    lon = coord[1]
    x = r * np.sin(col) * np.cos(lon)
    y = r * np.sin(col) * np.sin(lon)
    z = r * np.cos(col)
    return np.array([x, y, z])


def sphere2fish(sphere_point, r_pic):

    theta, phi = sphere_point
    r_f = r_pic * np.tan(theta/2)
    # r_f = r_pic * np.sin(theta)
    # r_f = 2 * r_pic * theta / np.pi
    # r_f = np.sqrt(2) * r_pic * np.sin(theta/2)
    x_f = r_f * np.cos(phi) + r_pic
    y_f = r_f * np.sin(phi) + r_pic
    return np.array([x_f, y_f])
