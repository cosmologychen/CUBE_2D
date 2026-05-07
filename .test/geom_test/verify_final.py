import numpy as np
import struct

# 使用修复后的读取逻辑 (= 前缀)
def load_lc_comparison(fn):
    import struct
    header_format = '=qfqqqqq' 
    header_size = struct.calcsize(header_format)
    with open(fn, 'rb') as f:
        header_data = f.read(header_size)
        res = struct.unpack(header_format, header_data)
        header = {'match_count': res[0], 'box': res[1], 'ng': res[2]}
        dt = np.dtype([('pid', np.int64), ('dist', np.float32), ('r_proj', np.float32), 
                       ('t_proj', np.float32), ('theta', np.float32)])
        data = np.fromfile(f, dtype=dt, count=header['match_count'])
    return header, data

header, data = load_lc_comparison('.test/geom_test/final_test.bin')
print(f"Match Count: {header['match_count']}")
# 检查前 5 个，看 dist 是否全为正，r_proj 是否有正有负
for i in range(5):
    print(f"PID {data['pid'][i]}: dist={data['dist'][i]:.6f}, r_proj={data['r_proj'][i]:.6f}")

# 验证 dist 是否有负数
neg_dist = np.any(data['dist'] < 0)
print(f"Has negative dist: {neg_dist}")
