import numpy as np
import struct

def load_fixed(fn):
    import struct
    header_format = '=qfqqqqq' 
    header_size = struct.calcsize(header_format)
    with open(fn, 'rb') as f:
        header_data = f.read(header_size)
        res = struct.unpack(header_format, header_data)
        mc = res[0]
        header = {'match_count': mc, 'box': res[1], 'ng': res[2]}
        data = np.zeros(mc, dtype=[('pid', 'i8'), ('dist', 'f4'), ('r_proj', 'f4'), ('t_proj', 'f4'), ('theta', 'f4')])
        data['pid'] = np.fromfile(f, dtype='i8', count=mc)
        data['dist'] = np.fromfile(f, dtype='f4', count=mc)
        data['r_proj'] = np.fromfile(f, dtype='f4', count=mc)
        data['t_proj'] = np.fromfile(f, dtype='f4', count=mc)
        data['theta'] = np.fromfile(f, dtype='f4', count=mc)
    return header, data

# 读取原始输入 A/B 计算真值
ng = 1024
pos_a = np.fromfile('.test/geom_test/A_fixed_xp.bin', dtype='f4').reshape(-1, 2)
pid_a = np.fromfile('.test/geom_test/A_fixed_pid.bin', dtype='i8')
pos_b = np.fromfile('.test/geom_test/B_fixed_xp.bin', dtype='f4').reshape(-1, 2)
pid_b = np.fromfile('.test/geom_test/B_fixed_pid.bin', dtype='i8')

map_a = {p: x for p, x in zip(pid_a, pos_a)}
map_b = {p: x for p, x in zip(pid_b, pos_b)}

header, data = load_fixed('.test/geom_test/final_validation.bin')
print(f"Validated {header['match_count']} matches.")

residuals = []
for i in range(header['match_count']):
    p = data['pid'][i]
    d_calc = data['dist'][i]
    dx = (map_b[p] - map_a[p] + ng/2.0) % ng - ng/2.0
    d_true = np.sqrt(np.sum(dx**2))
    residuals.append(abs(d_calc - d_true))

print(f"Max Residual: {np.max(residuals):.6e}")
# 检查 r_proj 是否有正有负
print(f"r_proj range: {np.min(data['r_proj']):.4f} to {np.max(data['r_proj']):.4f}")
# 检查 dist 是否全为正
print(f"dist min: {np.min(data['dist']):.4f}")
