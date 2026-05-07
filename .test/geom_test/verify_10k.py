import numpy as np
import struct

ng = 1024
pos_a = np.fromfile('.test/geom_test/A_10k_xp.bin', dtype=np.float32).reshape(-1, 2)
pid_a = np.fromfile('.test/geom_test/A_10k_pid.bin', dtype=np.int64)
pos_b = np.fromfile('.test/geom_test/B_10k_xp.bin', dtype=np.float32).reshape(-1, 2)
pid_b = np.fromfile('.test/geom_test/B_10k_pid.bin', dtype=np.int64)

map_a = {pid: pos for pid, pos in zip(pid_a, pos_a)}
map_b = {pid: pos for pid, pos in zip(pid_b, pos_b)}

with open('.test/geom_test/results_10k.bin', 'rb') as f:
    match_count = struct.unpack('q', f.read(8))[0]
    f.read(4 + 8 + 8*4) # 跳过头部其他信息
    
    residuals = []
    record_struct = struct.Struct('qffff')
    for _ in range(match_count):
        r = record_struct.unpack(f.read(record_struct.size))
        pid, d_calc = r[0], r[1]
        
        # 重新计算真值
        dx = (map_b[pid] - map_a[pid] + ng/2.0) % ng - ng/2.0
        d_true = np.sqrt(np.sum(dx**2))
        residuals.append(abs(d_calc - d_true))

print(f"Verified {match_count} matches.")
print(f"Max Residual: {np.max(residuals):.6e}")
print(f"Mean Residual: {np.mean(residuals):.6e}")
