import numpy as np

ng = 1024
np_base = 1000

# --- 生成列表 A ---
# 随机生成 1000 个粒子
pid_a = np.arange(1, np_base + 1, dtype=np.int64)
pos_a = np.random.uniform(0, ng, (np_base, 2)).astype(np.float32)

# --- 生成列表 B ---
# 1. 复制 A 的前 900 个粒子，并加入微小随机位移 (0.1 格点)
pid_b = pid_a[:900].copy()
pos_b = (pos_a[:900] + np.random.normal(0, 0.1, (900, 2))).astype(np.float32) % ng

# 2. 插入 50 个 B 独有的新粒子 (PID 从 2000 开始)
pid_b = np.concatenate([pid_b, np.arange(2000, 2050, dtype=np.int64)])
pos_b = np.concatenate([pos_b, np.random.uniform(0, ng, (50, 2)).astype(np.float32)])

# 3. 故意制造 10 个重复粒子 (取 A 中的前 10 个)
pid_b = np.concatenate([pid_b, pid_a[:10]])
pos_b = np.concatenate([pos_b, pos_a[:10]])

# 保存数据
pos_a.tofile('.test/geom_test/A_large_xp.bin')
pid_a.tofile('.test/geom_test/A_large_pid.bin')
pos_b.tofile('.test/geom_test/B_large_xp.bin')
pid_b.tofile('.test/geom_test/B_large_pid.bin')

print(f"Large test data generated: A={len(pid_a)}, B={len(pid_b)} (Expected: 900 match, 100 missing in B, 50 exclusive in B, 10 dups)")
