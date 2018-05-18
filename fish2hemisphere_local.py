# coding: utf-8
""" sphere_conv用のfish2hemisphere.py

"""

from __future__ import division, absolute_import, print_function

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image

SAVE_DIR = "./RESULT/spherical_array/"
if not os.path.isdir(SAVE_DIR):
    os.makedirs(SAVE_DIR)


def cal_value(channel, position, r_pic, im_in, pos_2d):

    channel = channel
    i_2d, j_2d = position
    x_f, y_f = sphere2fish(pos_2d[i_2d][j_2d], r_pic)
    # im_out[i, j, channel] = bilinear_interpolate(im_in[:, :, channel], x, y)
    value = bilinear_interpolate(im_in[:, :, channel], x_f, y_f)

    return value


def cal_value_from3d(channel, position, r_pic, im_in, pos_3d):
    """球面配列から，その位置における画素値をバイリニア補間で直接計算する．
    Parameters
    ----------
    channel : R, G, B チャンネルのこと (R:0, G:1, B:2)
    position : 球面配列における座標．(area, i, j)
    r_pic : 魚眼画像の画像部分の半径
    im_in : 魚眼画像のnparray
    pos_3d : 球面配列 Numpy配列orリスト (area, i, j, 2)みたいな配列

    Return
    ---------
    指定されたchannelにおける画素値
    """
    channel = channel
    area, i_3d, j_3d = position
    x_f, y_f = sphere2fish(pos_3d[area][i_3d][j_3d], r_pic)
    # im_out[i, j, channel] = bilinear_interpolate(im_in[:, :, channel], x, y)
    value = bilinear_interpolate(im_in[:, :, channel], x_f, y_f)

    return value

