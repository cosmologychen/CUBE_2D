import numpy as np
import struct

# --- 1. 读取原始生成的数据 A 和 B ---
# A 是基准点，B 是带扰动的点
pos_a_raw = np.fromfile('.test/geom_test/A_fixed_xp.bin', dtype=np.float32).reshape(-1, 2)
pid_a_raw = np.fromfile('.test/geom_test/A_fixed_pid.bin', dtype=np.int64)
pos_b_raw = np.fromfile('.test/geom_test/B_fixed_xp.bin', dtype=np.float32).reshape(-1, 2)
pid_b_raw = np.fromfile('.test/geom_test/B_fixed_pid.bin', dtype=np.int64)

# 建立 PID 到位置的映射，以便直接对比
map_a = {pid: pos for pid, pos in zip(pid_a_raw, pos_a_raw)}
map_b = {pid: pos for pid, pos in zip(pid_b_raw, pos_b_raw)}

# --- 2. 读取 Fortran 计算结果 ---
ng = 1024
with open('.test/geom_test/fixed_results.bin', 'rb') as f:
    def read_i8(): return struct.unpack('q', f.read(8))[0]
    def read_f4(): return struct.unpack('f', f.read(4))[0]

    match_count = read_i8()
    box = read_f4()
    ng_file = read_i8()
    f.read(8 * 4) # n_a, n_b, mis_a, mis_b

    record_struct = struct.Struct('qffff')
    calc_data = []
    for _ in range(match_count):
        calc_data.append(record_struct.unpack(f.read(record_struct.size)))

# --- 3. 直接计算真值距离并对比 ---
print(f"{'PID':<6} | {'True Dist':<10} | {'Calc Dist':<10} | {'Residual':<10}")
print("-" * 50)

residuals = []
for pid, dist_calc, r_p, t_p, theta in calc_data:
    pa = map_a[pid]
    pb = map_b[pid]
    
    # Python 计算 PBC 距离 (与 Fortran modulo 逻辑一致)
    dx = pb - pa
    dx = (dx + ng/2.0) % ng - ng/2.0
    d_true = np.sqrt(np.sum(dx**2))
    
    res = abs(dist_calc - d_true)
    residuals.append(res)
    
    if pid <= 10:
        print(f"{pid:<6} | {d_true:<10.6f} | {dist_calc:<10.6f} | {res:<10.6e}")

print("-" * 50)
print(f"平均残差 (Mean Residual): {np.mean(residuals):.6e}")
print(f"最大残差 (Max Residual):  {np.max(residuals):.6e}")
