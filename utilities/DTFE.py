import sys
import numpy as np
from scipy.spatial import Delaunay


def load2dxpos(fn):
    fid = open(fn, 'rb')
    p1 = np.fromfile(fid, dtype=np.float32)
    fid.close()
    n = round(len(p1)/2)
    if (n*2 != len(p1)):
        print('shape no mach')
        print(n,n**2,len(p1))
        return 0,0
    a = np.reshape(p1, (2, n),order='F')
    print(fn,n)
    return a.T

# 获取命令行参数
if len(sys.argv) < 2:
    print("Usage: python dtfe_extract.py <xp_filename>")
    sys.exit(1)

xp_file = sys.argv[1]

# 读取粒子数据
xp = load2dxpos(xp_file)

# Delaunay 三角剖分
tri = Delaunay(xp)
triangle_vertices = xp[tri.simplices]  # shape = (M, 3, 2)
# print(triangle_vertices.shape)
# print(np.max(triangle_vertices),np.min(triangle_vertices))

# 写入输出文件
triangle_vertices.astype(np.float32).tofile("triangles.bin")