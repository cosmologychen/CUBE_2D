import numpy as np
import struct

with open('.test/geom_test/test_results.bin', 'rb') as f:
    def read_i8(): return struct.unpack('q', f.read(8))[0]
    def read_f4(): return struct.unpack('f', f.read(4))[0]

    match_count = read_i8()
    box = read_f4()
    # 尝试读取 ng，如果 A 中的 PID 正确，说明这里可能没有填充或者需要跳过填充
    # 我们根据 PID 的值来判断读取是否偏移
    ng = read_i8()
    n_a = read_i8()
    n_b = read_i8()
    mis_a = read_i8()
    mis_b = read_i8()
    
    print(f"Header: Match={match_count}, Box={box}, ng={ng}")
    print(f"Stats: n_a={n_a}, n_b={n_b}, mis_a={mis_a}, mis_b={mis_b}")
    
    record_struct = struct.Struct('qffff')
    while True:
        chunk = f.read(record_struct.size)
        if len(chunk) < record_struct.size: break
        r = record_struct.unpack(chunk)
        print(f"PID {r[0]}: dist={r[1]:.2f}, r_proj={r[2]:.2f}, t_proj={r[3]:.2f}, theta={r[4]:.2f}")
