import numpy as np
np.random.seed(42)
ng = 1024
np_base = 1000
pid_a = np.arange(1, np_base + 1, dtype=np.int64)
pos_a = np.random.uniform(0, ng, (np_base, 2)).astype(np.float32)
delta = np.random.normal(0, 0.1, (900, 2)).astype(np.float32)
pos_b = (pos_a[:900] + delta) % ng
pid_b = pid_a[:900]
pos_a.tofile('.test/geom_test/A_fixed_xp.bin')
pid_a.tofile('.test/geom_test/A_fixed_pid.bin')
pos_b.tofile('.test/geom_test/B_fixed_xp.bin')
pid_b.tofile('.test/geom_test/B_fixed_pid.bin')
