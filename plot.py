# -*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D

def scatter3d(filename, sphere=True, view=(), nb_point=None, r=1):
    """3次元プロット
    Parameters
    ----------
    filename : 保存するfilename
    sphere : Trueにすると，北半球のグリッドを表示する
    view : 球面を見る角度を決められる
    nb_point : 指定すれば，指定した点の近傍を色付けして出力
    """
    ele, azi = view
    figure = plt.figure(figsize=(8, 8))
    ax = Axes3D(figure)
    # ax = figure.add_subplot(111, projection='3d')
    ax.set_aspect("equal")

    if sphere is True:
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = r * np.cos(u) * np.sin(v)
        y = r * np.sin(u) * np.sin(v)
        z = r * np.cos(v)
        ax.plot_wireframe(x, y, z, color='r', linewidth=0.5)

    cartesian_vertex = np.array([
        [spherical2cartesian(c, r=r) for c in point_vertex]
    ])
        cartesian_coords = np.array([
        [spherical2cartesian(c, r=r) for c in points]
    ])

    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.set_zlabel("Z axis")
    ax.set_xlim(-1.1*r, 1.1*r)
    ax.set_ylim(-1.1*r, 1.1*r)
    ax.set_zlim(-1.1*r, 1.1*r)
    ax.scatter3D(
        xs=cartesian_coords.T[0],
        ys=cartesian_coords.T[1],
        zs=cartesian_coords.T[2],
        linewidth=1,
        color='k',
        )
    ax.scatter3D(
        xs=cartesian_vertex.T[0],
        ys=cartesian_vertex.T[1],
        zs=cartesian_vertex.T[2],
        linewidth=5,
        color='b',
        )
    if nb_point is not None:
        area, j_cen, k_cen = nb_point
        neighbor = find_neighbor(point=nb_point)
        cartesian_neighbor = np.array([
            [spherical2cartesian(c, r=r) for c in neighbor]
        ])
        cartesian_center = spherical2cartesian(pos[area][j_cen][k_cen], r=r)
        ax.scatter3D(
            xs=cartesian_center[0],
            ys=cartesian_center[1],
            zs=cartesian_center[2],
            linewidth=8,
            color='m',
            marker='h'
            )
        ax.scatter3D(
            xs=cartesian_neighbor.T[0],
            ys=cartesian_neighbor.T[1],
            zs=cartesian_neighbor.T[2],
            linewidth=8,
            color='g',
            marker='*',
            alpha=1
            )

    ax.view_init(elev=ele, azim=azi)
    plt.savefig(SAVE_DIR+filename)
    plt.show()
    plt.close()
    return
