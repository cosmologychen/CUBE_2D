import numpy as np
import struct

ng = 1024
pos_a_raw = np.fromfile('.test/geom_test/A_10k_xp.bin', dtype=np.float32).reshape(-1, 2)
pid_a_raw = np.fromfile('.test/geom_test/A_10k_pid.bin', dtype=np.int64)
pos_b_raw = np.fromfile('.test/geom_test/B_10k_xp.bin', dtype=np.float32).reshape(-1, 2)
pid_b_raw = np.fromfile('.test/geom_test/B_10k_pid.bin', dtype=np.int64)

# 核心逻辑：模拟 Fortran 的去重，只取第一次出现的 PID
map_a = {}
for p, x in zip(pid_a_raw, pos_a_raw):
    if p not in map_a: map_a[p] = x

map_b = {}
for p, x in zip(pid_b_raw, pos_b_raw):
    if p not in map_b: map_b[p] = x

with open('.test/geom_test/results_10k.bin', 'rb') as f:
    match_count = struct.unpack('q', f.read(8))[0]
    f.read(4 + 8 + 8*4)
    
    record_struct = struct.Struct('qffff')
    residuals = []
    for _ in range(match_count):
        r = record_struct.unpack(f.read(record_struct.size))
        pid, d_calc = r[0], r[1]
        
        # 计算基于“第一次出现”的真值
        dx = (map_b[pid] - map_a[pid] + ng/2.0) % ng - ng/2.0
        d_true = np.sqrt(np.sum(dx**2))
        residuals.append(abs(d_calc - d_true))

print(f"Verified {len(residuals)} matches at 10k scale.")
print(f"Max Residual: {np.max(residuals):.6e}")
print(f"Mean Residual: {np.mean(residuals):.6e}")
