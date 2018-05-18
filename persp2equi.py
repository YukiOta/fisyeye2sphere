# -*- coding: utf-8 -*-

from tqdm import tqdm
from PIL import Image
from sklearn.cluster import KMeans
from multiprocessing import Pool
from glob import glob

from keras.datasets import mnist
# from keras.utils import np_utils

import numpy as np
import sys
import os
import argparse
import fish2hemisphere_local as fish
import image_transfer as it


# load argument
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_path', default='../data/',
                    help='Path to data dirctory')
parser.add_argument('-s', '--save_path', default='../data/',
                    help='Path to save dirctory')
parser.add_argument('-m', '--mnist', action="store_true",
                    help='Flag to try demo with MNIST data')
parser.add_argument('-r', '--rotation', type=int, default='x',
                    choices=['x', 'z', 'y', 'xz', 'xy'],
                    help='Direction of rotation (default: x)')
args = parser.parse_args()
data_path, save_path = args.data_path, args.save_path
if not os.path.isdir(save_path):
    print("Create new save dirctory at " + save_path)
    os.makedirs(save_path)


def main():
    # 画像の読み込み
    if args.mnist:
        (x_train, y_train), (x_test, y_test) = mnist.load_data()
    else:
        x_data = []
        x_list = glob(os.path.join(data_path, "*.png"))
        for path in tqdm(ls):
            img = Image.open(path)
            img = img.resize((100, 100))
            img = np.array(img, dtype=np.float32)
            x_data.append(img)
        x_data = np.array(x_data)

    img_ori = np.zeros((224, 224, 3))
    level = 4
    img_fish, pos = fish.make_considered_picture(img=img_ori, level=level, return_array=1)
    pos_x = rotation_xyz(pos, alpha=np.pi/2, mode="x")
    print("level: ", level)
    print(img_fish.shape)



if __name__ == '__main__':
    main()
