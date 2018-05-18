
def make_2d_array(pos_3d, Q):
    """3次元配列を2次元配列に変換する
    Parameters
    ----------
    pos_3d : area, j, k で構成されているリスト
    Q : 分割数
    """

    R_L = 3*Q
    i_2d, j_2d = (3*Q, 5*Q)
    pos_2d = np.zeros((i_2d+1, j_2d+1, 2))
    pos_2d_coord = np.zeros((i_2d+1, j_2d+1, 3), dtype=np.int)

    for i in range(R_L+1):
        if i == 0:  # 1
            pos_2d[i][0] = pos_3d[1][0][1]
            pos_2d_coord[i][0] = (1, 0, 1)
        elif i <= Q:  # 2
            j_p, k_p = (i, 1)
            for j in range(5*Q):
                point = j % i
                if j <= i-1:
                    pos_2d[i][j] = pos_3d[1][j_p-point][k_p+point]
                elif j <= 2*i-1:
                    pos_2d[i][j] = pos_3d[2][j_p-point][k_p+point]
                elif j <= 3*i-1:
                    pos_2d[i][j] = pos_3d[3][j_p-point][k_p+point]
                elif j <= 4*i-1:
                    pos_2d[i][j] = pos_3d[4][j_p-point][k_p+point]
                elif j <= 5*i-1:
                    pos_2d[i][j] = pos_3d[5][j_p-point][k_p+point]
        elif i <= Q*2:  # 3
            i_n = i - Q
            j_p, k_p = (Q, 1+i_n)
            for j in range(5*Q):
                point = j % Q
                if j <= Q-1:
                    pos_2d[i][j] = pos_3d[1][j_p-point][k_p+point]
                    pos_2d_coord[i][j] = (1, j_p-point, k_p+point)
                elif j <= 2*Q-1:
                    pos_2d[i][j] = pos_3d[2][j_p-point][k_p+point]
                    pos_2d_coord[i][j] = (2, j_p-point, k_p+point)
                elif j <= 3*Q-1:
                    pos_2d[i][j] = pos_3d[3][j_p-point][k_p+point]
                    pos_2d_coord[i][j] = (3, j_p-point, k_p+point)
                elif j <= 4*Q-1:
                    pos_2d[i][j] = pos_3d[4][j_p-point][k_p+point]
                    pos_2d_coord[i][j] = (4, j_p-point, k_p+point)
                elif j <= 5*Q-1:
                    pos_2d[i][j] = pos_3d[5][j_p-point][k_p+point]
                    pos_2d_coord[i][j] = (5, j_p-point, k_p+point)

        elif i < R_L:
            i_ = R_L - i
            j_p, k_p = (Q, 2*Q-i_)
            for j in range(5*i_):
                point = j % i_
                if j <= i_-1:
                    pos_2d[i][j] = pos_3d[1][j_p-point][k_p+point]
                    pos_2d_coord[i][j] = (1, j_p-point, k_p+point)
                elif j <= 2*i_-1:
                    pos_2d[i][j] = pos_3d[2][j_p-point][k_p+point]
                    pos_2d_coord[i][j] = (2, j_p-point, k_p+point)
                elif j <= 3*i_-1:
                    pos_2d[i][j] = pos_3d[3][j_p-point][k_p+point]
                    pos_2d_coord[i][j] = (3, j_p-point, k_p+point)
                elif j <= 4*i_-1:
                    pos_2d[i][j] = pos_3d[4][j_p-point][k_p+point]
                    pos_2d_coord[i][j] = (4, j_p-point, k_p+point)
                elif j <= 5*i_-1:
                    pos_2d[i][j] = pos_3d[5][j_p-point][k_p+point]
                    pos_2d_coord[i][j] = (5, j_p-point, k_p+point)
        elif i == R_L:
            pos_2d[i][0] = pos_3d[1][Q][2*Q+1]
            pos_2d_coord[i][0] = (1, Q, 2*Q+1)

    return pos_2d, pos_2d_coord


def make_2d_array_interpolated(pos_3d, Q):
    """3次元配列を2次元配列に変換する (for シータ・ファイ画像)
    theta_phi画像を作る為に，天頂付近の座標は，細かく分割する．
    Parameters
    ----------
    pos_3d : area, j, k で構成されているリスト
    Q : 分割数
    """

    R_L = 3*Q
    i_2d, j_2d = (3*Q, 5*Q)
    pos_2d = np.zeros((i_2d, j_2d, 2))
    pos_2d_coord = np.zeros((i_2d, j_2d, 3), dtype=np.int)
    # phi_delta = (2*np.pi) / (5*Q)
    tmp = np.linspace(0, 2*np.pi, 5*Q, endpoint=True)

    for i in range(R_L+1):
        if i == 0:  # 1
            pos_2d[i] = pos_3d[1][0][1]
            pos_2d_coord[i] = (1, 0, 1)
        elif i <= Q:  # 2
            theta = pos[1][i][1][0]
            pos_2d[i, :, 0] = theta
            pos_2d[i, :, 1] = tmp
        elif i <= Q*2:  # 3
            i_n = i - Q
            j_p, k_p = (Q, 1+i_n)
            if pos[1][j_p][k_p][0] <= np.pi/2:
                theta = pos[1][j_p][k_p][0]
                pos_2d[i, :, 0] = theta
                pos_2d[i, :, 1] = tmp

    return pos_2d


def make_2d_picture(pos_2d, r_pic, Q, im_in):

    cols, rows = (pos_2d.shape[1], pos_2d.shape[0])
    im_out = np.zeros((rows, cols, 3), np.uint8)

    R_L = 3*Q
    for channel in range(3):
        for i in range(R_L+1):
            if i == 0:  # 1
                for j in range(5):
                    im_out[i, j, channel] = cal_value(channel, position=(i, 0), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i <= Q:  # 2
                for j in range(5*i):
                    im_out[i, j, channel] = cal_value(channel, position=(i, j), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i <= Q*2:  # 3
                if pos_2d[i][0][0] <= np.pi/2:
                    for j in range(5*Q):
                        im_out[i, j, channel] = cal_value(channel, position=(i, j), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i < R_L:
                if pos_2d[i][0][0] <= np.pi/2:
                    i_ = R_L - i
                    for j in range(5*i_):
                        im_out[i, j, channel] = cal_value(channel, position=(i, j), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i == R_L:
                if pos_2d[i][0][0] <= np.pi/2:
                    im_out[i, 0, channel] = cal_value(channel, position=(i, 0), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)

    return im_out[:int(im_out.shape[0]/2), :, :]


def make_2d_theta_phi_picture(pos_2d, r_pic, Q, im_in):

    cols, rows = (pos_2d.shape[1], pos_2d.shape[0])
    im_out = np.zeros((rows, cols, 3), np.uint8)

    R_L = 3*Q
    for channel in range(3):
        for i in range(rows):
            for j in range(cols):
                im_out[i, j, channel] = cal_value(channel, position=(i, j), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)

    return im_out[:int(im_out.shape[0]/2), :, :]


def make_2d_center_picture(pos_2d, r_pic, Q, im_in):

    cols, rows = (pos_2d.shape[1], pos_2d.shape[0])
    im_out = np.zeros((rows, cols, 3), np.uint8)

    R_L = 3*Q
    for channel in range(3):
        for i in range(R_L+1):
            if i == 0:  # 1
                index = int(5*Q / 2 - 1)
                im_out[i, index, channel] = cal_value(channel, position=(i, 0), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i <= Q:  # 2
                if i % 2 == 0:
                    index = int((5*Q - 5*i) / 2)
                else:
                    index = int((5*Q - 5*i - 1) / 2)
                for j in range(5*i):
                    im_out[i, index+j, channel] = cal_value(channel, position=(i, j), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i <= Q*2:  # 3
                if pos_2d[i][0][0] <= np.pi/2:
                    for j in range(5*Q):
                        im_out[i, j, channel] = cal_value(channel, position=(i, j), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i < R_L:
                if pos_2d[i][0][0] <= np.pi/2:
                    i_ = R_L - i
                    for j in range(5*i_):
                        im_out[i, j, channel] = cal_value(channel, position=(i, j), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i == R_L:
                if pos_2d[i][0][0] <= np.pi/2:
                    im_out[i, 0, channel] = cal_value(channel, position=(i, 0), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)

    return im_out


def make_2d_round_picture(pos_2d, r_pic, Q, im_in):
    """画像に途切れる部分が内容に，画素値を，行において繰り返す関数
    Parameters
    ----------
    """
    cols, rows = (pos_2d.shape[1], pos_2d.shape[0])
    im_out = np.zeros((rows, cols-1, 3), np.uint8)

    R_L = 3*Q
    for channel in range(3):
        for i in range(R_L+1):
            value_tmp = []
            if i == 0:  # 1
                im_out[i, :, channel] = cal_value(channel, position=(i, 0), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i <= Q:  # 2
                if i % 2 == 0:
                    index = int((5*Q - 5*i) / 2)
                else:
                    index = int((5*Q - 5*i - 1) / 2)
                for j in range(5*i):
                    tmp = cal_value(channel, position=(i, j), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
                    value_tmp.append(tmp)
                for k in range(5*i+index+1):
                    if index+k < 80:
                        im_out[i, index+k, channel] = value_tmp[k % (5*i)]
                    if index-k >= 0:
                        im_out[i, index-k, channel] = value_tmp[4-((k-1) % (5*i))]
            elif i <= Q*2:  # 3
                if pos_2d[i][0][0] <= np.pi/2:
                    for j in range(5*Q):
                        im_out[i, j, channel] = cal_value(channel, position=(i, j), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i < R_L:
                if pos_2d[i][0][0] <= np.pi/2:
                    i_ = R_L - i
                    for j in range(5*i_):
                        im_out[i, j, channel] = cal_value(channel, position=(i, j), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i == R_L:
                if pos_2d[i][0][0] <= np.pi/2:
                    im_out[i, 0, channel] = cal_value(channel, position=(i, 0), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)

    return im_out[:int(im_out.shape[0]/2), :, :]


def make_2d_neighbor_considered_picture(pos_2d, pos_2d_coord, r_pic, Q, im_in):
    """画像に途切れる部分が内容に，find_neighborを用いて，一つ下の列から
    画素値を補間
    Parameters
    ----------
    pos_2d_coord : 2次元配列の座標に対する，(area, i, j)が格納されてる
    neib_coord[1] : 対象とするpointの右上の点
    """
    cols, rows = (pos_2d.shape[1], pos_2d.shape[0])
    im_out = np.zeros((rows, cols, 3), np.uint8)

    for channel in range(3):
        for i in reversed(range(Q*2+1)):
            if i == 0:  # 1
                im_out[i, :, channel] = cal_value(channel, position=(i, 0), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
            elif i <= Q:  # 2
                for j in range(5*Q):
                    neib_pos, neib_coord = find_neighbor(point=pos_2d_coord[i+1, j], coord=True, level=level)
                    area_, i_, j_ = neib_coord[1]
                    pos_2d_coord[i, j] = neib_coord[1]
                    im_out[i, j, channel] = img_sphere[channel, area_, i_, j_]
            elif i <= Q*2:  # 3
                if pos_2d[i][0][0] <= np.pi/2:
                    for j in range(5*Q):
                        im_out[i, j, channel] = cal_value(channel, position=(i, j), r_pic=r_pic, im_in=im_in, pos_2d=pos_2d)
    return im_out[:int(im_out.shape[0]/2), :, :]


def make_considered_picture(img, level, return_array=0):
    """Parameters
    r : 分割レベル
    Q : 平行四辺形の短辺における点の数．座標にもなる．
    i : 平行四辺形 (area) の数 (に，実装の都合上+1してる)
    j : 平行四辺形の短辺の座標の上限 (に，実装の都合上+1してる)
    k : 平行四辺形の長辺の座標の上限 (に，実装の都合上+1してる)
    kai : bc方向(円周方向)の単位角度 (論文参照)
    tau : ab方向(天頂から天底方向)の単位角度
    img_sphere : 球面配列に画素値が入るような配列 (3, i, j, k)
    """
    global r_pic, im_in, r, Q, i, j, k, kai, tau, img_sphere, pos, points
    global point_vertex
    r = level
    Q = 2**r
    i = 5 + 1
    j = Q + 1 + 1
    k = 2*Q + 1 + 1
    kai = 2 * np.pi / 5
    tau = np.arctan(2)
    img_sphere = np.zeros((3, i, j, k))

    """画像の読み込み
    PILを使ってImage.open
    ここでは計算の不可を少なくするため，(100, 100)にリサイズしている

    r_pic : 魚眼画像の半径 (画像がある部分)
    """
    # im_in = Image.open("../1data/test_image/cloud_half.png")
    # im_in = im_in.resize((100, 100))
    # im_in = np.array(im_in)
    im_in = img
    r_pic = int(im_in.shape[0]/2)


    """球面座標を保持するリストの定義
    pos[area][短辺座標][長辺座標]

    Parameters
    ----------
    pos[i][j][k][0] : ido theta [0, pi]
    pos[i][j][k][1] : keido phi [0, 2*pi]
    """
    pos = np.zeros((i, j, k, 2))
    pos[1][0][1] = (0, 0)  # a
    pos[1][Q][1] = (tau, 0)  # b
    pos[1][0][Q+1] = (tau, kai)  # c
    pos[1][Q][Q+1] = (np.pi - tau, kai / 2)  # d
    pos[1][0][2*Q+1] = (np.pi - tau, kai * 3 / 2)  # e
    pos[1][Q][2*Q+1] = (np.pi, 0)  # f

    """点のリスト
    points : 分割された点を放り込んでいく
    point_vertex : 頂点の座標を放り込んでいく
    vertex_list : area1の球面座標は手打ちでいれている (初期値的なね)
    """
    points = []
    point_vertex = []
    vertex_list = [
        (1, 0, 1),
        (1, Q, 1),
        (1, 0, Q+1),
        (1, Q, Q+1),
        (1, 0, 2*Q+1),
        (1, Q, 2*Q+1),
    ]

    """北半球にある点群だけを取り出してくる
    # vertex : こっちでは，area1を計算してる
    # ALL vertexes of icosahedron : こっちでは，全ての面を計算する
    """
    # vertex
    for point in vertex_list:
        area, j_v, k_v = point
        check_hemisphere(
            pos=pos,
            point=(area, j_v, k_v), p_list=point_vertex,
            img_sphere=img_sphere, make_img_array=True
        )

    # ALL vertexes of icosahedron
    for point in vertex_list:
        area, j_v, k_v = point
        for i in range(1, 5):
            pos[i+1][j_v][k_v] = cal_angle(pos[i][j_v][k_v])
            check_hemisphere(
                pos=pos,
                point=(i+1, j_v, k_v), p_list=point_vertex,
                img_sphere=img_sphere, make_img_array=True
            )

    """平行四辺形の分割
    基本的に，areaの値以外変えていない
    Parameters
    ----------
    area : どの平行四辺形か
    a-f : 頂点の座標を指定する

    returnは無いが，暗示的に先ほど定義したpointsリストに座標を放り込んでいる
    """
    for i in range(1, 6):
        tessellation(area=i, a=(0, 1), b=(Q, 1), c=(0, Q+1), d=(Q, Q+1), e=(0, 2*Q+1), f=(Q, 2*Q+1))
        # print(len(points))

    """これは実験的な，area1のみの分割に対応する
    """
    # tessellation(area=1, a=(0, 1), b=(Q, 1), c=(0, Q+1), d=(Q, Q+1), e=(0, 2*Q+1), f=(Q, 2*Q+1))

    """点の可視化
    Parameters
    ----------
    filename : 保存するfilenameを指定できる
    view=(a, b) : どの角度で見るか
    area : area
    nb_point=(area, i, j) : 指定すれば，指定した点の周りの近傍を色付けして出力する
    r : 魚眼画像の円部分の半径である．
    """
    a, b, area, i, j = (30, 30, 1, Q, 1)
    # scatter3d(filename='hemi_%d_%d_%d_%d%d%d.png' % (r, a, b, area, i, j), view=(a, b), nb_point=(area, i, j), r=r_pic)
    # scatter3d(filename='hemi_%d_%d_%d_%d%d%d.png' % (r, a, b, area, i, j), view=(a, b), nb_point=None, r=r_pic)

    """球面上に配置した点の座標から，球面画像の作成
    """
    pos_2d, pos_2d_coord = make_2d_array(pos_3d=pos, Q=Q)
    # pos_2d = make_2d_array_interpolated(pos_3d=pos, Q=Q)
    if return_array == 1:
        return img_sphere, pos
    # im = make_2d_picture(pos_2d=pos_2d, r_pic=r_pic, Q=Q, im_in=im_in)
    # im_center = make_2d_center_picture(pos_2d=pos_2d, r_pic=r_pic, Q=Q, im_in=im_in)
    # im_consider = make_2d_neighbor_considered_picture(pos_2d_coord=pos_2d_coord, pos_2d=pos_2d, r_pic=r_pic, Q=Q, im_in=im_in, level=level)

    return im_consider

