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
    f.read(4 + 8 + 8*4)
    
    worst_res = -1.0
    worst_data = None
    
    record_struct = struct.Struct('qffff')
    for _ in range(match_count):
        r = record_struct.unpack(f.read(record_struct.size))
        pid, d_calc = r[0], r[1]
        
        pa = map_a[pid]
        pb = map_b[pid]
        
        # 尝试两种 Python PBC 逻辑
        # 1. 直接用 %
        dx_py = (pb - pa + ng/2.0) % ng - ng/2.0
        d_true = np.sqrt(np.sum(dx_py**2))
        
        res = abs(d_calc - d_true)
        if res > worst_res:
            worst_res = res
            worst_data = (pid, pa, pb, d_calc, d_true, dx_py)

print(f"Worst Residual: {worst_res:.6e}")
pid, pa, pb, dc, dt, dxp = worst_data
print(f"PID {pid}:")
print(f"  Pos A: {pa}")
print(f"  Pos B: {pb}")
print(f"  Calc: {dc:.6f}, True: {dt:.6f}")
print(f"  Python DX: {dxp}")
