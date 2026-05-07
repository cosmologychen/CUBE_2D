import numpy as np

pid_b = np.fromfile('.test/geom_test/B_10k_pid.bin', dtype=np.int64)
pos_b = np.fromfile('.test/geom_test/B_10k_xp.bin', dtype=np.float32).reshape(-1, 2)

# 找出 PID 60 出现的所有位置
indices = np.where(pid_b == 60)[0]
for idx in indices:
    print(f"Index {idx}, PID 60 Pos: {pos_b[idx]}")

