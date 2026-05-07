# CUBE_2D 程序结构文档

## 1. 概述

CUBE_2D 是一个 2D N-body  cosmological simulation 程序，使用 Coarray Fortran 并行化和 FFTW 进行快速傅里叶变换。程序模拟宇宙中暗物质粒子的演化，采用 PM (Particle-Mesh) 方法结合 PP (Particle-Particle) 校正计算引力。

## 2. 编译与运行

### 编译
```bash
make        # 编译主程序和工具
```

### 运行流程 (run.sh)
```bash
./utilities/ic.x      # 生成初始条件
./main.x               # 主模拟循环
./utilities/density_power.x  # 计算功率谱
./utilities/dsp.x      # 位移场分析
./utilities/fof.x      # Friends-of-Friends 星系团查找
```

## 3. 程序模块结构

### 3.1 核心模块

#### parameters.f90 - 参数模块
定义所有模拟参数：

| 参数 | 值 | 说明 |
|------|-----|------|
| ndim | 2 | 维度 |
| box | 600 | 盒子大小 (Mpc/h) |
| ng | 3072 | 网格分辨率 |
| ncore | 12 | OpenMP 线程数 |
| ratio_cs | 4 | 粗细网格比 |
| nc | ng/ratio_cs = 768 | 粗网格分辨率 |
| nnt | 4 | 每图像每维瓦片数 |
| nt | nc/nnt = 192 | 每瓦片粗网格数 |
| np_nc | ratio_cs | 每粗网格粒子数 |
| apm1c | 3.5 | PM1 软化长度 |
| apm2 | 3.5 | PM2 软化长度 |
| app | 0.06 | PP 软化长度 |
| omega_m | 0.309 | 物质密度参数 |
| omega_l | 0.691 | 暗能量密度参数 |
| istep_max | 100000 | 最大时间步数 |

**关键类型定义:**
- `sim_header`: 存储模拟元信息 (np, a, t, dt, checkpoint等)
- `type_halo_catalog_header/array`: 星系团目录结构

#### variables.f90 - 变量模块
声明全局变量和辅助子程序：

**粒子数据:**
- `xp(2, np_max)`: 粒子位置
- `vp(2, np_max)`: 粒子速度
- `pid(np_max)`: 粒子ID（在 2D 模拟中，`pid(ip) = ip`，主要用于兼容 3D 架构）

**密度数组:**
- `rho1(nw+2, nw)`: 细网格密度 (与rho1k等价)
- `rho1k(nw/2+1, nw)`: 细网格密度 (傅里叶空间)
- `rho2(1-ngb:ngp+ngb, 1-ngb:ngp+ngb, ncore)`: 粗网格密度
- `rho2k(ngt/2+1, ngt, ncore)`: 粗网格密度 (傅里叶空间)

**FFT计划:**
- `plan, iplan`: 细网格FFT计划
- `plan2(ncore), iplan2(ncore)`: 粗网格FFT计划

**关键子程序:**
- `spine_tile()`: 计算瓦片上粒子累积索引
- `spine_image()`: 计算图像上粒子累积索引

#### Green.f90 - 格林函数
`Green_2D()`: 计算2D格林函数用于PM计算

#### basic_functions.f08 - 基础函数
- `tic/toc`: 计时函数
- `wrap_position()`: 周期性边界条件处理
- `pbc_vec()`: 周期性向量差
- `xpos2mesh()`: 位置转网格索引
- `output_name()`: 生成输出文件名
- `F_ra()`: PP力计算中的径向函数

### 3.2 主程序模块

#### initialize.f90
初始化程序：
1. 设置OpenMP线程
2. 读取 `z_checkpoint.txt` 红移列表
3. 计算/读取格林函数 Gk1, Gk2
4. 创建FFT计划
5. 初始化PP邻居列表

#### particle_initialization.f90
读取检查点数据：
1. 读取 `info` 文件获取sim结构
2. 读取 `xp`, `vp` 文件获取粒子位置和速度

#### timestep.f90
时间步进计算：
- `expansion()`: 根据宇宙学模型计算 da/dt
- 自适应调整dt，确保膨胀率 ra < ra_max
- 检查是否需要checkpoint

#### drift.f90
位置更新（drift操作）：
```fortran
xp = xp + dt_mid * vp
```

#### kick.f90
速度更新（kick操作），包含三种力计算：

**PM1 (粗网格引力):**
- 将粒子云密度分配到细网格 (CIC)
- FFT计算势能
- 插值力到粒子位置

**PM2 (细网格引力):**
- 分瓦片(tile)处理，每个瓦片有缓冲区
- FFT计算势能
- 更新粒子速度

**PP (粒子对相互作用):**
- 使用链表 (hoc/ll) 组织粒子
- 对邻居单元格中的粒子计算PP力
- 软化长度: app = 0.06

**时间步长约束:**
```fortran
dt = min(dt_e, dt_pm1, dt_pm2, dt_refine*dt_pp, dt_vmax)
```

#### checkpoint.f90
保存检查点：
- 写入 `info`: sim结构
- 写入 `xp`: 粒子位置
- 写入 `vp`: 粒子速度
- 写入 `pid`: 粒子ID (索引)

#### runtime_lightcone.f90
RunTime 光锥检测模块：
- **预测性检测**：在 `drift` 之前利用 `dt_mid * vp` 预判粒子下一位置。
- **插值逻辑**：支持一阶线性插值和二阶 Newton-Raphson 迭代。
- **并行写入**：使用线程局部缓冲区 (`thread_buffers`) 避免并行锁。

#### main.f90
主循环：
```fortran
do istep = timestep, istep_max
    call timestep      ! 计算dt
    
    ! 光锥预测检测
    if (enable_runtime_lightcone) then
        call check_lightcone_crossing(xp, vp, pid, sim%np)
    endif
    
    call drift         ! 更新位置
    call kick          ! 更新速度
    if (checkpoint_step) then
        call checkpoint
    endif
enddo
```

包含宇宙学生长因子函数 `Dgrow()`。

## 4. 工具程序 (utilities/)

#### ic.f90 - 初始条件生成
1. 读取功率谱数据
2. 生成高斯随机噪声
3. Wiener滤波
4. 计算势能并获取位移场
5. 生成粒子初始位置和速度
6. 输出: `xp`, `vp`, `info`, `delta_L`, `phi1`

**输出文件:**
- `delta_L.bin`: 线性密度场
- `phi1.bin`: 初始势能
- `xp.bin`, `vp.bin`: 粒子位置速度
- `info.bin`: 模拟信息

#### density_power.f90 - 密度功率谱
1. CIC密度分配到网格
2. FFT计算功率谱
3. 互相关功率谱计算 (可选)

**输出:** `delta_c.bin`, `power.bin`

#### displacement.f90 - 位移场分析
计算和分析 Lagrangian 位移场：
- `dsp_D`: 原始位移场
- `dsp_sD`: 平滑后位移场
- `dsp_E`: 发散位移场
- `delta_E`: 发散密度
- 分解位移场为膨胀/剪切/旋转分量

#### displacement_sm.f90 - 平滑位移场分析
使用高斯平滑的位移场分析版本。

#### displacement_DTFE.f90 - DTFE位移场分析
使用 DTFE (Delaunay Tessellation Field Estimator) 方法进行位移场分析。

#### fof.f90 - Friends-of-Friends 星系团查找
Friends-of-Friends算法：
- 连接长度 b_link = 0.20
- 最小粒子数 np_halo_min = 2
- 使用链表组织网格中的粒子

**输出:** `halo.bin`, `halo_xp_mean_only.bin`, `halo_qp_mean_only.bin`

#### void.f90 - 空洞查找
基于低密度区域的空洞 finder。

#### DT_void.f90 - DTFE空洞查找
使用 Delaunay Tessellation 方法查找空洞。

#### powerspectrum.f90 - 功率谱计算模块
包含功率谱计算的相关子程序：
- `cross_power()`: 互功率谱
- `auto_power()`: 自功率谱
- `pk_correction()`: 功率谱校正

#### test_lc_geom.f90 - 光锥几何比对工具
比对不同算法生成的光锥粒子：
1. **PID 匹配**：基于原始索引进行精确配对。
2. **空间哈希**：利用粗网格链表 (HOC/LL) 实现 $O(N)$ 匹配。
3. **几何分解**：计算位移偏差、径向投影 $R$ 和切线投影 $\Theta_R$。

## 5. 数据文件格式

### 检查点文件
| 文件 | 内容 | 格式 |
|------|------|------|
| `{z}_info.bin` | sim_header结构 | stream |
| `{z}_xp.bin` | 粒子位置 | stream, real(4) |
| `{z}_vp.bin` | 粒子速度 | stream, real(4) |
| `{z}_pid.bin` | 粒子ID (索引) | stream, int(8) |

### 光锥输出文件 (RunTime / Interp)
| 文件 | 内容 |
|------|------|
| `*_xp.bin` | 光锥面上粒子的插值位置 |
| `*_vp.bin` | 穿过瞬间的粒子速度 |
| `*_pid.bin` | 原始粒子索引 (用于比对) |
| `*_a.bin` | 穿过瞬间的尺度因子 |

### 分析输出文件
| 文件 | 内容 |
|------|------|
| `{z}_delta_c.bin` | 密度反差 |
| `{z}_power.bin` | 功率谱 |
| `{z}_dsp_D.bin` | 位移场 |
| `{z}_halo.bin` | 星系团目录 |
| `lc_comparison_results.txt` | 光锥几何比对结果 |

## 6. 关键算法

### 6.1 PM算法 (Particle-Mesh)
1. **密度分配**: 使用CIC (Cloud-in-Cell) 将粒子分配到网格
2. **FFT**: 快速傅里叶变换到k空间
3. **格林函数**: 与格林函数卷积计算势能
4. **力计算**: 有限差分法从势能计算力
5. **速度更新**: 插值力到粒子位置并更新速度

### 6.2 多层网格结构
```
粗网格 (nc=768) ──> PM1 (范围: ratio_cs*apm1c = 14)
细网格 (ng=3072) ──> PM2 (范围: apm2 = 3.5)
PP校正 (范围: app = 0.06)
```

### 6.3 周期性边界条件
使用 `wrap_position()` 和 `pbc_vec()` 处理周期性边界。

### 6.4 链表组织 (hoc/ll)
- `ll(ip)`: 粒子ip的下一个粒子索引
- 支持快速邻居查找

### 6.5 光锥跨越检测与插值 (Lightcone Crossing)
1. **预测性判定**：计算 $r_i$ 和 $r_{i+1}$（预测），当 $(r_i - \chi_i)(r_{i+1} - \chi_{i+1}) \le 0$ 时判定跨越。
2. **插值比例**：计算穿越比例 $f = d_i / (d_i - d_{i+1})$。
3. **状态插值**：
   - 一阶：$x_c = x_i + f(x_{i+1} - x_i)$。
   - 二阶：利用 $c_{eff} = -a(d\chi/dt)$ 进行 Newton-Raphson 修正。

## 7. 并行化

- **OpenMP**: 主要并行化方式
- **Coarray Fortran**: 用于多图像通信（本版本未完全使用）
- **FFTW**: 多线程FFT计算

## 8. 宇宙学参数

| 参数 | 值 |
|------|-----|
| h0 | 0.6766 |
| omega_m | 0.309 |
| omega_l | 0.691 |
| omega_nu | 0 ( neutrinos) |
| n_s | 0.9665 |
| sigma_8 | 0.821 |

## 9. 修改程序注意事项

### 添加新参数
1. 在 `parameters.f90` 的适当位置添加
2. 如需在运行时读取，添加到 `sim_header` 类型
3. 在 `initialize.f90` 或 `checkpoint.f90` 中处理读写

### 添加新模块变量
1. 在 `variables.f90` 中声明
2. 注意内存使用，大数组考虑分配

### 修改引力计算
- PM1: `kick.f90` 的 "PM1" 部分
- PM2: `kick.f90` 的 "PM2" 部分
- PP: `kick.f90` 的 "PP" 部分

### 添加新工具程序
1. 在 `utilities/` 目录创建
2. 修改 `utilities/Makefile` 添加编译规则
3. 包含 `use parameters` 和必要的模块

### 输出修改
- 使用 `output_name()` 函数生成带红移的文件名
- 文件保存在 `opath` 目录
