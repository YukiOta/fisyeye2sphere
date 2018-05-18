# -*- coding: utf-8 -*-

def cal_angle(pos):
    """角度計算する関数
    正二十面体の各頂点を計算するときに使ってる
    """
    theta = pos[0]
    # phi = (kai + pos[1]) % (2 * np.pi)
    phi = (kai + pos[1])
    return (theta, phi)


def cal_mid(area=1, coord1=(0, 0), coord2=(1, 1)):
    """中点を計算するプログラム
    Parameters
    ----------

    """

    j_1, k_1 = coord1
    j_2, k_2 = coord2
    j = int((j_1 + j_2) / 2)
    k = int((k_1 + k_2) / 2)

    if coord1 == (0, 1):
        mid_0 = (pos[area][j_1][k_1][0] + pos[area][j_2][k_2][0]) / 2
        mid_1 = pos[area][j_2][k_2][1]
        mid_point = (mid_0, mid_1)
    elif coord2 == (0, 1):
        mid_0 = (pos[area][j_1][k_1][0] + pos[area][j_2][k_2][0]) / 2
        mid_1 = pos[area][j_1][k_1][1]
        mid_point = (mid_0, mid_1)
    else:
        mid_point = (pos[area][j_1][k_1] + pos[area][j_2][k_2]) / 2
    return (j, k, mid_point)


def cal_unit_12(area=1, coord1=(0, 0), coord2=(1, 1)):
    """1 -> 2の方向への単位ベクトルを計算する関数
    vertexの間の単位ベクトルを計算するときに使っているっぽい？
    Parameters
    ----------
    """
    j_1, k_1 = coord1
    j_2, k_2 = coord2
    # j = int((j_1 + j_2) / 2)
    # k = int((k_1 + k_2) / 2)

    if coord1 == (0, 1):
        start = pos[area][j_1][k_1][0]
        end = pos[area][j_2][k_2][0]
        edge_unit = np.array(((end - start) / Q, 0))
    elif coord2 == (Q, 2*Q+1):
        start = pos[area][j_1][k_1][0]
        end = pos[area][j_2][k_2][0]
        edge_unit = np.array(((end - start) / Q, 0))
    else:
        start = pos[area][j_1][k_1]
        end = pos[area][j_2][k_2]
        edge_unit = (end - start) / Q

    return edge_unit


def cal_midunit_12(area=1, coord1=(0, 0), coord2=(1, 1), sep=2):
    """真ん中の部分での単位ベクトル？？
    cal_unit_12との違いは，パラメータにsepを含む子こと．
    使用時は，for文でsepの値を指定しているから，こっちを使っている
    Parameters
    ----------
    area : どの平行四辺形か (1-5)
    coord1 : 始点
    coord2 : 終点
    sep : 何分割するかのパラメータ
    """
    j_1, k_1 = coord1
    j_2, k_2 = coord2
    # j = int((j_1 + j_2) / 2)
    # k = int((k_1 + k_2) / 2)

    if coord1 == (0, 1):
        start = pos[area][j_1][k_1][0]
        end = pos[area][j_2][k_2][0]
        edge_unit = np.array(((end - start) / sep, 0))
    elif coord2 == (Q, 2*Q+1):
        start = pos[area][j_1][k_1][0]
        end = pos[area][j_2][k_2][0]
        edge_unit = np.array(((end - start) / sep, 0))
    else:
        start = pos[area][j_1][k_1]
        end = pos[area][j_2][k_2]
        edge_unit = (end - start) / sep

    return edge_unit


def tessellation(
        area=1, a=(0, 1), b=(4, 1), c=(0, 5),
        d=(0, 5), e=(4, 5), f=(0, 9)):
    """球面を分割する関数
    Parameters
    -----------
    area : 平行四辺形の番号 (1~5)
    a - f : 平行四辺形を最初に三角形で分割したときの，拡張点
    """
    j_a, k_a = a
    j_b, k_b = b
    j_c, k_c = c
    j_d, k_d = d
    j_e, k_e = e
    j_f, k_f = f

    unit_ab = cal_unit_12(area=area, coord1=(j_a, k_a), coord2=(j_b, k_b))
    unit_ac = cal_unit_12(area=area, coord1=(j_a, k_a), coord2=(j_c, k_c))
    unit_bc = cal_unit_12(area=area, coord1=(j_b, k_b), coord2=(j_c, k_c))
    unit_cd = cal_unit_12(area=area, coord1=(j_c, k_c), coord2=(j_d, k_d))
    unit_ce = cal_unit_12(area=area, coord1=(j_c, k_c), coord2=(j_e, k_e))
    unit_de = cal_unit_12(area=area, coord1=(j_d, k_d), coord2=(j_e, k_e))
    unit_bd = cal_unit_12(area=area, coord1=(j_b, k_b), coord2=(j_d, k_d))
    unit_df = cal_unit_12(area=area, coord1=(j_d, k_d), coord2=(j_f, k_f))
    unit_ef = cal_unit_12(area=area, coord1=(j_e, k_e), coord2=(j_f, k_f))

    for i in range(1, Q):
        # ABC
        pos[area][i][k_a] = np.array(
            (pos[area][j_a][k_a][0], pos[area][j_b][k_b][1])) + i * unit_ab
        pos[area][j_a][k_a+i] = np.array(
            (pos[area][j_a][k_a][0], pos[area][j_c][k_c][1])) + i * unit_ac
        pos[area][j_b-i][k_b+i] = pos[area][j_b][k_b] + i * unit_bc
        # CDE
        pos[area][i][k_c] = pos[area][j_c][k_c] + i * unit_cd
        pos[area][j_c][k_c+i] = pos[area][j_c][k_c] + i * unit_ce
        pos[area][j_d-i][k_d+i] = pos[area][j_d][k_d] + i * unit_de
        # BD
        pos[area][j_b][k_b+i] = pos[area][j_b][k_b] + i * unit_bd
        # DEF
        pos[area][j_d][k_d+i] = pos[area][j_d][k_d] + i * unit_df
        pos[area][j_e+i][k_e] = pos[area][j_e][k_e] + i * unit_ef

        # check_hemisphere(pos=pos[area][i][k_a], p_list=points)
        # check_hemisphere(pos=pos[area][j_a][k_a+i], p_list=points)
        # check_hemisphere(pos=pos[area][j_b-i][k_b+i], p_list=points)
        #
        # check_hemisphere(pos=pos[area][i][k_c], p_list=points)
        # check_hemisphere(pos=pos[area][j_c][k_c+i], p_list=points)
        # check_hemisphere(pos=pos[area][j_d-i][k_d+i], p_list=points)
        #
        # check_hemisphere(pos=pos[area][j_b][k_b+i], p_list=points)
        # check_hemisphere(pos=pos[area][j_d][k_d+i], p_list=points)
        # check_hemisphere(pos=pos[area][j_e+i][k_e], p_list=points)
        check_hemisphere(
            pos=pos, point=(area, i, k_a), p_list=points,
            img_sphere=img_sphere, make_img_array=True
        )
        check_hemisphere(
            pos=pos, point=(area, j_a, k_a+i), p_list=points,
            img_sphere=img_sphere, make_img_array=True
        )

        check_hemisphere(
            pos=pos, point=(area, j_b-i, k_b+i), p_list=points,
            img_sphere=img_sphere, make_img_array=True
        )

        check_hemisphere(
            pos=pos, point=(area, i, k_c), p_list=points,
            img_sphere=img_sphere, make_img_array=True
        )
        check_hemisphere(
            pos=pos, point=(area, j_c, k_c+i), p_list=points,
            img_sphere=img_sphere, make_img_array=True
        )
        check_hemisphere(
            pos=pos, point=(area, j_d-i, k_d+i), p_list=points,
            img_sphere=img_sphere, make_img_array=True
        )

        check_hemisphere(
            pos=pos, point=(area, j_b, k_b+i), p_list=points,
            img_sphere=img_sphere, make_img_array=True
        )
        check_hemisphere(
            pos=pos, point=(area, j_d, k_d+i), p_list=points,
            img_sphere=img_sphere, make_img_array=True
        )
        check_hemisphere(
            pos=pos, point=(area, j_e+i, k_e), p_list=points,
            img_sphere=img_sphere, make_img_array=True
        )

    # make points points inside of triangles
    # from left to right, horizontal direction
    for i in range(1, Q):
        if i >= 2:
            unit_abc = cal_midunit_12(area=area, coord1=(i, k_a), coord2=(j_a, k_a+i), sep=i)
            unit_cde = cal_midunit_12(area=area, coord1=(i, k_c), coord2=(j_c, k_c+i), sep=i)
            for j in range(1, i):
                pos[area][i-j][k_a+j] = pos[area][i][k_a] + j * unit_abc
                pos[area][i-j][k_c+j] = pos[area][i][k_c] + j * unit_cde

                check_hemisphere(
                    pos=pos, point=(area, i-j, k_a+j), p_list=points,
                    img_sphere=img_sphere, make_img_array=True
                )
                check_hemisphere(
                    pos=pos, point=(area, i-j, k_c+j), p_list=points,
                    img_sphere=img_sphere, make_img_array=True
                )

                # check_hemisphere(pos=pos[area][i-j][k_a+j], p_list=points)
                # check_hemisphere(pos=pos[area][i-j][k_c+j], p_list=points)

        unit_bcd = cal_midunit_12(area=area, coord1=(j_b, k_b+i), coord2=(j_c+i, k_c), sep=Q-i)
        unit_def = cal_midunit_12(area=area, coord1=(j_d, k_d+i), coord2=(j_e+i, k_e), sep=Q-i)
        for j in range(1, Q-i):
            pos[area][j_b-j][k_b+i+j] = pos[area][j_b][k_b+i] + j * unit_bcd
            pos[area][j_d-j][k_d+i+j] = pos[area][j_d][k_d+i] + j * unit_def
            check_hemisphere(
                pos=pos, point=(area, j_b-j, k_b+i+j), p_list=points,
                img_sphere=img_sphere, make_img_array=True
            )
            check_hemisphere(
                pos=pos, point=(area, j_d-j, k_d+i+j), p_list=points,
                img_sphere=img_sphere, make_img_array=True
            )

            # check_hemisphere(pos=pos[area][j_b-j][k_b+i+j], p_list=points)
            # check_hemisphere(pos=pos[area][j_d-j][k_d+i+j], p_list=points)

    return


def check_hemisphere(pos, point, p_list, img_sphere, make_img_array=True):
    """半球に抑えこむためのチェッカー
    pos[0]には緯度の情報が入ってる

    一緒に，画像も作ってしまう

    Parameters
    ----------
    pos : 球面配列が入ったリスト
    point : (area, j, k)
    p_list : 北半球にあれば，p_listにその座標を追加する
    img_sphere : 最初の方で定義しておく，画素情報を入れ込む配列
    """
    area, j__, k__ = point
    channels = 3
    '''check sphereにしちゃう
    '''
    if pos[area][j__][k__][0] <= np.pi:
        p_list.append(pos[area][j__][k__])
        if make_img_array is True:
            for channel in range(channels):
                img_sphere[channel][area, j__, k__] = cal_value_from3d(
                    channel=channel,
                    position=(area, j__, k__),
                    r_pic=r_pic,
                    im_in=im_in,
                    pos_3d=pos
                )
        return


def find_neighbor(level, point=(1, 1, 2), coord=False):
    """注目点の近傍を持ってくる
    Parameters
    ----------
    point : 注目する点 (area, j ,k)

    return
    ----------
    neighbor : neighborの曲座標
    neighbor_coord : neigborの球面座標
    """
    Q = 2**level
    # print('in find_neighbor')
    # print(Q)
    area, j_p, k_p = check_neighbor(level=level, point=point)
    neighbor = []
    neighbor_coord = []
    neighbor_list = [
        (area, j_p-1, k_p+1),
        (area, j_p-1, k_p),
        (area, j_p, k_p-1),
        (area, j_p+1, k_p-1),
        (area, j_p+1, k_p),
        (area, j_p, k_p+1)
    ]
    if (j_p, k_p) == (0, 1):  # zenith
        for i in range(1, 6):
            neighbor.append(pos[i][1][1])
            neighbor_coord.append((i, 1, 1))
        '''zenithのとき、真ん中の値を足す
        '''
        neighbor.append(pos[1][0][1])
        neighbor_coord.append((1, 0, 1))
    elif (j_p, k_p) == (Q, 2*Q+1):  # nadir
        for i in range(1, 6):
            neighbor.append(pos[i][Q][2*Q])
            neighbor_coord.append((i, Q, 2*Q))
        '''nadirのときも、真ん中値を足す
        '''
        neighbor.append(pos[1][Q][2*Q+1])
        neighbor_coord.append((1, j_p, k_p))
    else:
        for point in neighbor_list:
            area, j_nb, k_nb = point
            area, j_nb, k_nb = check_neighbor(level=level, point=(area, j_nb, k_nb))
            neighbor.append(pos[area][j_nb][k_nb])
            neighbor_coord.append((area, j_nb, k_nb))

    if coord is True:
        return neighbor, neighbor_coord
    else:
        return neighbor


def check_neighbor(level, point=(1, 2, 1)):
    """近傍のチェッカー
    Parameters
    ----------
    point : 注目する点 (area, j, k)
    """
    Q = 2**level
    # print('in check_neighbor')
    # print(Q)
    area, j, k = point
    if k == 0:
        area = area - 1
        if area == 0:
            area = 5
        return (area, 1, j)

    elif j == Q+1:
        area = area - 1
        if area == 0:
            area = 5
        if k <= Q:
            return (area, 1, Q+k)
        else:
            return (area, k-Q, 2*Q)

    elif j == 0:
        area = area + 1
        if area == 6:
            area = 1
        if k >= 2 and k <= Q:
            return (area, k-1, 1)
        elif k >= Q+1 and k <= 2*Q+1:
            return (area, Q, k-Q)
        else:
            return (area, j, k)

    elif k == 2*Q+1:
        area = area + 1
        if area == 6:
            area = 1
        return (area, Q, j+Q+1)

    else:
        return (area, j, k)
