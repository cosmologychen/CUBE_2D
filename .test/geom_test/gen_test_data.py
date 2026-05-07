import numpy as np

ng = 1024

# List A
pos_a = np.array([[100.0, 100.0], [500.0, 500.0], [10.0, 10.0], [1023.5, 512.0]], dtype=np.float32)
pid_a = np.array([1, 2, 3, 4], dtype=np.int64)

# List B
pos_b = np.array([[100.0, 100.0], [501.0, 500.0], [10.0, 10.0], [10.0, 10.0], [0.5, 512.0], [200.0, 200.0]], dtype=np.float32)
pid_b = np.array([1, 2, 3, 3, 4, 5], dtype=np.int64)

pos_a.tofile('.test/geom_test/A_xp.bin')
pid_a.tofile('.test/geom_test/A_pid.bin')
pos_b.tofile('.test/geom_test/B_xp.bin')
pid_b.tofile('.test/geom_test/B_pid.bin')

print("Synthetic test data generated in .test/geom_test/")
