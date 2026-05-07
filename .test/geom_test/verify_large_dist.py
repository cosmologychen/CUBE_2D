import numpy as np
import struct

# --- 1. 读取真值 (生成测试数据时加入的扰动) ---
# 在之前的脚本中，前 900 个粒子位移是 np.random.normal(0, 0.1, (900, 2))
# 我们需要重新生成并固定这个扰动以便比对
np.random.seed(42) # 固定种子
ng = 1024
np_match = 900
pos_a = np.random.uniform(0, ng, (1000, 2)).astype(np.float32)
# 真值位移 (注意处理 PBC)
delta_true = np.random.normal(0, 0.1, (np_match, 2)).astype(np.float32)
dist_true = np.sqrt(np.sum(delta_true**2, axis=1))

# --- 2. 读取 Fortran 计算结果 ---
with open('.test/geom_test/large_results.bin', 'rb') as f:
    def read_i8(): return struct.unpack('q', f.read(8))[0]
    def read_f4(): return struct.unpack('f', f.read(4))[0]

    match_count = read_i8()
    box = read_f4()
    ng_file = read_i8()
    f.read(8 * 4) # 跳过 n_a, n_b, mis_a, mis_b

    record_struct = struct.Struct('qffff')
    calc_results = []
    for _ in range(match_count):
        chunk = f.read(record_struct.size)
        calc_results.append(record_struct.unpack(chunk))

# --- 3. 比对 A 列表前 10 个匹配项的残差 ---
print(f"{'PID':<6} | {'True Dist':<10} | {'Calc Dist':<10} | {'Residual':<10}")
print("-" * 45)

residuals = []
for i in range(match_count):
    pid, dist_calc, r_proj, t_proj, theta = calc_results[i]
    if pid <= np_match: # 对应我们设置的有扰动的粒子
        d_true = dist_true[pid-1]
        res = abs(dist_calc - d_true)
        residuals.append(res)
        if i < 10: # 仅打印前 10 个示例
            print(f"{pid:<6} | {d_true:<10.6f} | {dist_calc:<10.6f} | {res:<10.6e}")

print("-" * 45)
print(f"平均残差 (Mean Residual): {np.mean(residuals):.6e}")
print(f"最大残差 (Max Residual):  {np.max(residuals):.6e}")
