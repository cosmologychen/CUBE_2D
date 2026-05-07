import numpy as np

ng = 1024
np_base = 10000
np_match = 9500

# --- A 列表 ---
pid_a = np.arange(1, np_base + 1, dtype=np.int64)
pos_a = np.random.uniform(0, ng, (np_base, 2)).astype(np.float32)

# --- B 列表 ---
# 1. 匹配部分 (9500个)
pid_b = pid_a[:np_match].copy()
pos_b = (pos_a[:np_match] + np.random.normal(0, 0.2, (np_match, 2))).astype(np.float32) % ng

# 2. B 独有部分 (1000个)
pid_b = np.concatenate([pid_b, np.arange(50000, 51000, dtype=np.int64)])
pos_b = np.concatenate([pos_b, np.random.uniform(0, ng, (1000, 2)).astype(np.float32)])

# 3. 重复部分 (100个)
pid_b = np.concatenate([pid_b, pid_a[:100]])
pos_b = np.concatenate([pos_b, pos_a[:100]])

pos_a.tofile('.test/geom_test/A_10k_xp.bin')
pid_a.tofile('.test/geom_test/A_10k_pid.bin')
pos_b.tofile('.test/geom_test/B_10k_xp.bin')
pid_b.tofile('.test/geom_test/B_10k_pid.bin')
print(f"10k Data Ready: A={len(pid_a)}, B={len(pid_b)}")
