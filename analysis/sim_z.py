# %%
import camb, os, struct, re, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
import matplotlib as mpl

# %% [markdown]
# ### 1. 基础配置与工具函数

# %%

plt.rcParams.update({
    'figure.figsize': (12, 8),
    'figure.facecolor': 'w',
    'figure.dpi': 300,
    'lines.linewidth': 1.0,
    'axes.spines.top': True,
    'axes.spines.right': True,
    'xtick.bottom': True,
    'xtick.direction': 'in',
    'ytick.left': True,
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'axes.linewidth': 1.0,
    'axes.xmargin': 0.03,
    'axes.ymargin': 0.03,
    'axes.grid': False,
    'axes.grid.which': 'both',
    'axes.labelpad': 8.0,
    'axes.labelsize': 10,
    'axes.labelcolor': 'k',
    'axes.axisbelow': True,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'xtick.minor.size': 2,
    'ytick.minor.size': 2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'xtick.minor.width': 1.0,
    'ytick.minor.width': 1.0,
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'font.size': 10,
    'text.usetex': False,
    'legend.fontsize': 12,
    'legend.frameon': False,
    'lines.markersize': 6,
    'axes.prop_cycle': plt.cycler('color', [[0, 0, 0], [1, 0.2, 0.2], [0.4, 0.7, 1], [0.2, 0.8, 0.3], [1, 0.7, 0.3], [0.6, 0.3, 1], [1, 0.7, 0.7], [0.5, 0.5, 0.5]])
})

def match_para(para):
    file_path = '../parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()
    pattern = r'parameter\s*::\s*%s\s*=\s*([^\s!]+)'%para
    match = re.search(pattern, content)
    if match:
        variable_value = match.group(1).strip().replace("'","").replace('"','')
        try:
            return float(variable_value)
        except:
            return variable_value
    else:
        print(para,"Pattern not found in the file.")
        sys.exit()

def get_sim_info(prefix):
    sim = {}
    with open(prefix + 'info.bin', 'rb') as fid:
        sim['np'] = struct.unpack('q', fid.read(8))[0]
        sim['izipx'] = struct.unpack('q', fid.read(8))[0]
        sim['izipv'] = struct.unpack('q', fid.read(8))[0]
        sim['nnt'] = struct.unpack('q', fid.read(8))[0]
        sim['nt'] = struct.unpack('q', fid.read(8))[0]
        sim['ncell'] = struct.unpack('q', fid.read(8))[0]
        sim['ncb'] = struct.unpack('q', fid.read(8))[0]
        sim['istep'] = struct.unpack('q', fid.read(8))[0]
        sim['cur_checkpoint'] = struct.unpack('q', fid.read(8))[0]
        sim['cur_powerpoint'] = struct.unpack('q', fid.read(8))[0]
        sim['calculate_PK'] = struct.unpack('q', fid.read(8))[0]
        sim['cic_iapm'] = struct.unpack('q', fid.read(8))[0]
        sim['a'] = struct.unpack('f', fid.read(4))[0]
        sim['t'] = struct.unpack('f', fid.read(4))[0]
        sim['tau'] = struct.unpack('f', fid.read(4))[0]
        sim['dt'] = struct.unpack('4f', fid.read(16))
        sim['mass_p'] = struct.unpack('f', fid.read(4))[0]
        sim['m_nu'] = struct.unpack('3f', fid.read(12))
        sim['Mass_nu'] = struct.unpack('f', fid.read(4))[0]
        sim['box'] = struct.unpack('f', fid.read(4))[0]
        sim['h0'] = struct.unpack('f', fid.read(4))[0]
        sim['omega_m'] = struct.unpack('f', fid.read(4))[0]
        sim['omega_l'] = struct.unpack('f', fid.read(4))[0]
        sim['s8'] = struct.unpack('f', fid.read(4))[0]
        sim['vsim2phys'] = struct.unpack('f', fid.read(4))[0]
        sim['sigma_vres'] = struct.unpack('f', fid.read(4))[0]
        sim['sigma_vi'] = struct.unpack('f', fid.read(4))[0]
        sim['z_i'] = struct.unpack('f', fid.read(4))[0]
        sim['vz_max'] = struct.unpack('f', fid.read(4))[0]
    sim['nc'] = sim['nt'] * sim['nnt']
    sim['ng_global'] = sim['nc'] * 4 
    return sim

def loadfield2d(fn):
    fid = open(fn, 'rb')
    p1 = np.fromfile(fid, dtype=np.float32)
    fid.close()
    n = round(len(p1) ** (1/2))
    a = np.reshape(p1, (n, n), order='F')
    return a, n

def load2dxpos(fn):
    fid = open(fn, 'rb')
    p1 = np.fromfile(fid, dtype=np.float32)
    fid.close()
    n = round(len(p1)/2)
    a = np.reshape(p1, (2, n), order='F')
    return a, int(np.sqrt(n))

def loadpower(filename):
    n_row_xi = 10
    fid = open(filename, 'rb')
    xi = np.fromfile(fid, dtype='float32')
    fid.close()
    xi = np.reshape(xi, (int(len(xi) / n_row_xi), n_row_xi))
    ksim = xi[:, 1]
    return ksim, xi

def load_lc_comparison(fn):
    """
    读取 test_lc_geom.x 生成的二进制比对结果
    返回: header (dict), data (structured array)
    """
    import struct
    # 头部格式: match_count(q=i8), box(f=f4), ng, n_a, n_b, mis_a, mis_b (5*q=i8)
    # 关键修复: 添加 '=' 前缀禁用 C 语言的自动对齐填充，确保与 Fortran 的 stream 严格对应 52 字节
    header_format = '=qfqqqqq' 
    header_size = struct.calcsize(header_format)
    
    if not os.path.exists(fn):
        print(f"File not found: {fn}")
        return None, None
        
    with open(fn, 'rb') as f:
        header_data = f.read(header_size)
        res = struct.unpack(header_format, header_data)
        header = {
            'match_count': res[0],
            'box': res[1],
            'ng': res[2],
            'n_a': res[3],
            'n_b': res[4],
            'mis_a': res[5],
            'mis_b': res[6]
        }
        # 记录格式: 由于 Fortran 采用批量数组写入，数据在文件中按字段块分布
        # 顺序: match_count个PID, 然后是match_count个dist, r_proj, t_proj, theta
        dt = np.dtype([('pid', np.int64), ('dist', np.float32), ('r_proj', np.float32), 
                       ('t_proj', np.float32), ('theta', np.float32)])

        data = np.zeros(header['match_count'], dtype=dt)
        data['pid']    = np.fromfile(f, dtype=np.int64,   count=header['match_count'])
        data['dist']   = np.fromfile(f, dtype=np.float32, count=header['match_count'])
        data['r_proj'] = np.fromfile(f, dtype=np.float32, count=header['match_count'])
        data['t_proj'] = np.fromfile(f, dtype=np.float32, count=header['match_count'])
        data['theta']  = np.fromfile(f, dtype=np.float32, count=header['match_count'])

        return header, data

# %%
# 初始化路径与参数
Path = '/mnt/18T/output_2D/lc_c10/3000_1024_iclc/'
z_list = np.loadtxt('../z_checkpoint.txt')
Redshift = f"{z_list[0]:.3f}"
sim = get_sim_info(Path + Redshift + '_')
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

# %% [markdown]
# ### 1.5 宇宙膨胀历史 (Expansion History: z vs chi)

# %%
eh_file = os.path.join(Path, 'neutrinos/s_a_tau_H.txt')
if os.path.exists(eh_file):
    eh_data = np.loadtxt(eh_file)
    # 根据 Pk.py 写入顺序: t, a_ex, tau, chi_ex
    eh_t, eh_a, eh_tau, eh_chi = eh_data
    eh_z = 1.0 / eh_a - 1.0
    
    # 仅保留 z 在 0-2 范围内的数据，以保证 Y 轴缩放正确
    mask = (eh_z >= 0) & (eh_z <= 2.0)
    plot_z = eh_z[mask]
    plot_chi = np.abs(eh_chi[mask])
    
    plt.figure(figsize=(8, 5), dpi=300)
    plt.plot(plot_z, plot_chi, label='Comoving Distance $\chi(z)$', color='red', linewidth=2)
    plt.xlabel('Redshift $z$')
    plt.ylabel('$\chi$ (Internal Units)')
    plt.title('Expansion History: $z$ vs $\chi$ (z=0-2)')
    plt.grid(True, alpha=0.3)
    # 标注模拟红移检查点 (仅显示 z <= 2)
    for checkpoint_z in z_list:
        if checkpoint_z <= 2.0:
            plt.axvline(x=checkpoint_z, color='gray', linestyle=':', alpha=0.5)
    plt.xlim(0, 2)
    plt.legend()
    plt.show()
    
    print(f"Max Redshift in EH: {eh_z.max():.4f}")
    print(f"Max Chi at max z: {np.abs(eh_chi).max():.4e}")
else:
    print(f"Expansion history file not found at {eh_file}")


# %%
eh_file = os.path.join(Path, 'neutrinos/s_a_tau_H.txt')
if os.path.exists(eh_file):
    eh_data = np.loadtxt(eh_file)
    # 根据 Pk.py 写入顺序: t, a_ex, tau, chi_ex
    eh_t, eh_a, eh_tau, eh_chi = eh_data
    eh_z = 1.0 / eh_a - 1.0
    
    # 仅保留 z 在 0-2 范围内的数据，以保证 Y 轴缩放正确
    mask = (eh_z >= 0) & (eh_z <= 2.0)
    plot_t = eh_t[mask]
    plot_chi = np.abs(eh_chi[mask])
    
    plt.figure(figsize=(8, 5), dpi=300)
    plt.plot(plot_t, plot_chi, label='Comoving Distance $\chi(z)$', color='red', linewidth=2)
    plt.xlabel('Time $t$')
    plt.ylabel('$\chi$ (Internal Units)')
    plt.grid(True, alpha=0.3)
    # 标注模拟红移检查点 (仅显示 t <= 2)
    for checkpoint_z in z_list:
        if checkpoint_z <= 2.0:
            plt.axvline(x=checkpoint_z, color='gray', linestyle=':', alpha=0.5)
    # plt.xlim(0, 2)
    plt.legend()
    plt.show()
    
    print(f"Max Time in EH: {eh_t.max():.4f}")
    print(f"Max Chi at max t: {np.abs(plot_chi).max():.4e}")
else:
    print(f"Expansion history file not found at {eh_file}")

# %% [markdown]
# ### 2. plt ic (1x4 子图)

# %%
end = 10
title = 'box = %d'%match_para('box')
pos, n_side = load2dxpos(f'{Path}{Redshift}_xp.bin')
k, xi = loadpower(f'{Path}' + Redshift + '_power.bin')
[k_camb, pk_camb] = [np.loadtxt(f'{Path}neutrinos/IC/Pcb_ic.txt')[:,0], np.loadtxt(f'{Path}neutrinos/IC/Pcb_ic.txt')[:,1]]

plt.figure(figsize=(16, 4), dpi=300)

# # 子图 1: 粒子散点
# plt.subplot(1, 4, 1)
# plt.scatter(pos[0, ::100], pos[1, ::100], c='black', s=.1, marker='o', edgecolors='none', alpha=0.5)
# plt.axis('equal')
# plt.axis('off')
# plt.title('Subsampled xps')

# 子图 2: P(k) 对比
plt.subplot(1, 4, 2)
plt.loglog(k[:-end], xi[:-end, 2], label=r'$P(k)$', color=colors[2], linewidth=2)
plt.loglog(k_camb, pk_camb, '--', label=r'$Pk_{CAMB}$', color=colors[0], linewidth=2)
plt.xlim(np.nanmin(k[:-end])/1.2, np.nanmax(k[:-end])*1.2)
plt.title(title)
plt.grid(True, 'both')
plt.legend()

# 子图 3: delta_L 密度场
plt.subplot(1, 4, 3)
delta, n = loadfield2d(f'{Path}{Redshift}_delta_L.bin')
plt.imshow(delta, cmap='gray', origin='lower')
plt.title(r'$\delta_L$')
plt.colorbar()
plt.axis('equal')
plt.axis('off')

# 子图 4: phi_L 势场
plt.subplot(1, 4, 4)
delta, n = loadfield2d(f'{Path}{Redshift}_phi1.bin')
plt.imshow(delta, cmap='gray', origin='lower')
plt.title(r'$\phi_L$')
plt.colorbar()
plt.axis('equal')
plt.axis('off')

plt.tight_layout()
plt.show()

# %% [markdown]
# ### 3. xps (Particle Positions)

# %%
# plt xpos s
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
import matplotlib as mpl
mpl.rcParams['animation.embed_limit'] = 1000  # 单位是 MB
frame_ids = list(range(0, 30))
[f_l,f_h,f_d] = [4,3,300]
fig, ax = plt.subplots(figsize=(f_l,f_h),dpi=f_d)

ng = sim['ng_global']
placeholder = np.zeros((ng, ng))
sc = ax.imshow(placeholder,cmap='gray', origin='lower') 
ax.axis('equal')
ax.axis('off')
sc.set_clim(-1,3)
plt.colorbar(sc)
title = ax.set_title('')
ax.set_aspect('equal', adjustable='box')

def update(i):
    delta,n = loadfield2d(f'{Path}{z_list[i]:.3f}_delta_c.bin')
    sc.set_data(delta.T)
    title.set_text(f'z = {z_list[i]:.3f}')
    return sc, title

ani = FuncAnimation(fig, update, frames=frame_ids, interval=300)

HTML(ani.to_jshtml())  # 在Jupyter中显示动画

# %% [markdown]
# ### 4. 功率谱 (Power Spectrum Evolution)

# %%
# 绘制所有红移的功率谱演化
plt.figure(figsize=(10, 6), dpi=300)
for i, z in enumerate(z_list):
    k, xi = loadpower(f'{Path}{z:.3f}_power.bin')
    plt.loglog(k[:-end], xi[:-end, 2], label=f'$z={z:.3f}$', color=colors[i % len(colors)], linewidth=2)
plt.loglog(k_camb, pk_camb, 'k--', alpha=0.5, label='CAMB IC')
plt.xlim(np.nanmin(k[:-end])/1.2, np.nanmax(k[:-end])*1.2)
plt.grid(True, 'both')
plt.xlabel('k [h/Mpc]')
plt.ylabel('P(k)')
plt.title('Power Spectrum Evolution')
plt.legend(ncol=3, fontsize=8)
plt.show()

# 绘制生长率比值
plt.figure(figsize=(10, 6), dpi=300)
k0, xi0 = loadpower(f'{Path}{z_list[0]:.3f}_power.bin')
for i, z in enumerate(z_list[::2]):
    k, xi = loadpower(f'{Path}{z:.3f}_power.bin')
    growth_ratio = xi[:-end, 2] / xi0[:-end, 2] / (((z_list[0]+1)/(z+1))**2)
    plt.loglog(k[:-end], growth_ratio, label=f'$z={z:.3f}$', color=colors[i % len(colors)], linewidth=2)
plt.axhline(1.0, color='k', linestyle='--')
plt.title('Growth Factor Ratio')
plt.legend()
plt.show()

# %% [markdown]
# ### 5. 互功率谱 (Cross Power Spectrum)

# %%
plt.figure(figsize=(8, 6), dpi=300)
for i, z in enumerate(z_list[::2]):
    k, xiLN = loadpower(f'{Path}{z:.3f}_Cpower_LN.bin')
    r_LN = xiLN[:, 8] / np.sqrt(xiLN[:, 5] * xiLN[:, 2])
    plt.plot(k, r_LN, label=f'$z={z:.3f}$')
plt.axhline(y=0.6, color='k', linestyle='--')
plt.axvline(x=0.1, color='k', linestyle='--')
plt.axhline(y=0.5, color='b', linestyle='--')
plt.axvline(x=0.4, color='b', linestyle='--')
plt.axhline(y=0.8, color='g', linestyle='--')
plt.axvline(x=0.03, color='g', linestyle='--')
plt.xscale('log')
plt.ylim(-0.01, 1.1)
plt.xlabel(r'$k$ [$h$/Mpc]')
plt.ylabel(r'$r_{LN}$')
plt.title('Cross Correlation Coefficient')
plt.legend()
plt.show()

# %% [markdown]
# ### 6. 密度场 (Density Field Plots & Animation)

# %%
# 通用密度场绘图函数 (传入文件名)
def plot_field_from_file(fn, title='', clim=[-1, 3], cmap='gray'):
    if not os.path.exists(fn):
        print(f"File not found: {fn}")
        return
    delta, n = loadfield2d(fn)
    plt.figure(figsize=(6, 5), dpi=300)
    plt.imshow(delta.T, cmap=cmap, origin='lower')
    plt.title(title)
    if clim is not None:
        plt.clim(clim[0], clim[1])
    plt.colorbar()
    plt.axis('off')
    plt.show()

# 模拟快照密度场绘图函数 (传入红移)
def plot_snapshot_field(z, name='delta_c', clim=[-1, 3], cmap='gray'):
    filename = f'{Path}{z:.3f}_{name}.bin'
    title = f'z = {z:.3f} ({name})'
    plot_field_from_file(filename, title, clim=clim, cmap=cmap)

# 调用展示示例
plot_snapshot_field(z_list[0])
plot_snapshot_field(z_list[-1])

# %%
# 动画展示
mpl.rcParams['animation.embed_limit'] = 1000 
fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
placeholder = np.zeros((sim['nc']*4, sim['nc']*4)) 
sc = ax.imshow(placeholder, cmap='gray', origin='lower') 
ax.axis('off')
sc.set_clim(-1, 3)
plt.colorbar(sc)
title = ax.set_title('')

def update(i):
    delta, n = loadfield2d(f'{Path}{z_list[i]:.3f}_delta_c.bin')
    sc.set_data(delta.T)
    title.set_text(f'z = {z_list[i]:.3f}')
    return sc, title

ani = FuncAnimation(fig, update, frames=range(len(z_list)), interval=300)
HTML(ani.to_jshtml())





# %% 切片光锥密度场检验
delta_c0 , n = loadfield2d(Path+'0.000_delta_c.bin')
delta_c06 , n = loadfield2d(Path+'0.687_delta_c.bin')
delta_lc , n = loadfield2d(Path+'lightcone_post_m1_delta.bin')

clim=[-1, 3]
cmap='gray'
[f_l,f_h,f_d] = [4,3,300]
plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow(delta_c0.T, cmap=cmap, origin='lower')
plt.title(f'z = 0.000 ({r"$\delta_c$"})')
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()


plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow(delta_c06.T, cmap=cmap, origin='lower')
plt.title(f'z = 0.687 ({r"$\delta_c$"})')
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()

plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow(delta_lc.T, cmap=cmap, origin='lower')
plt.title(f'z = 0.000 ({r"$\delta_{lc}$"})')
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()


plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow((delta_c0-delta_lc).T, cmap=cmap, origin='lower')
plt.title(f'z = 0.000 ({r"$\delta_c-\delta_{lc}$"})')
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()

plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow((delta_c06-delta_lc).T, cmap=cmap, origin='lower')
plt.title(f'z = 0.687 ({r"$\delta_c-\delta_{lc}$"})')
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()


# %% 光锥初始条件
Path_org = '/mnt/18T/output_2D/lc_c10/3000_1024_runtime/'
Path_lc = '/mnt/18T/output_2D/lc_c10/3000_1024_iclc/'
delta_L_org , n = loadfield2d(Path_org+'200.000_delta_L.bin')
delta_L_lc , n = loadfield2d(Path_lc+'200.000_delta_L.bin')
clim=[-1, 3]
cmap='gray'
[f_l,f_h,f_d] = [4,3,300]
plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow(delta_L_org.T, cmap=cmap, origin='lower')
plt.title(r"$\delta_{L,org}$")   
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()

plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow(delta_L_lc.T, cmap=cmap, origin='lower')
plt.title(r"$\delta_{L,lc}$")   
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()

plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow(delta_L_lc.T/delta_L_org.T, cmap=cmap, origin='lower')
plt.title(r"$\delta_{L,lc}-\delta_{L,org}$")   
# plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()

# %% 切片光锥密度场检验
delta_c0 , n = loadfield2d(Path_org+'0.000_delta_c.bin')
delta_c_lcic , n = loadfield2d(Path_lc+'0.000_delta_c.bin')
delta_lc , n = loadfield2d(Path_org+'lightcone_post_m1_delta.bin')  

clim=[-1, 3]
cmap='gray'
[f_l,f_h,f_d] = [4,3,300]
plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow(delta_c0.T, cmap=cmap, origin='lower')
plt.title(f'z = 0.000 ({r"$\delta_c$"})')
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()


plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow(delta_c_lcic.T, cmap=cmap, origin='lower')
plt.title(f'z = 0.000 ({r"$\delta_{c,lcic}$"})')
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()

plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow(delta_lc.T, cmap=cmap, origin='lower')
plt.title(f'z = 0.000 ({r"$\delta_{lc}$"})')
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()


plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow((delta_c0-delta_lc).T, cmap=cmap, origin='lower')
plt.title(f'z = 0.000 ({r"$\delta_c-\delta_{lc}$"})')
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()

plt.figure(figsize=(f_l,f_h),dpi=f_d)
plt.imshow((delta_c_lcic-delta_lc).T, cmap=cmap, origin='lower')
plt.title(f'z = 0.000 ({r"$\delta_{c,lcic}--\delta_{lc}$"})')
plt.clim(clim[0], clim[1])
plt.colorbar()
plt.axis('off')
plt.show()
# %%

# %% [markdown]
# ### 7. 光锥粒子比对可视化 (Lightcone Particle Comparison)

# %%
header, res = load_lc_comparison('/mnt/18T/output_2D/lc_c10/comparsion/runtime_runtime_2_runtime_slide.bin')
if res is not None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=200)
    fig.suptitle(f"Lightcone Comparison Statistics\n(Match: {header['match_count']}, Missing A: {header['mis_b']}, Missing B: {header['mis_a']}) \n Runtime .vs. Slide ", fontsize=14)

    # 1. 绝对距离分布
    axes[0, 0].hist(res['dist'], bins=100, color='skyblue', edgecolor='black', alpha=0.7)
    axes[0, 0].set_title(r'Absolute Distance Distribution (Grid Units)')
    axes[0, 0].set_xlabel(r'$\Delta r$')
    axes[0, 0].set_yscale('log')
    axes[0, 0].grid(True, alpha=0.3)

    # 2. 径向偏移分布
    axes[0, 1].hist(res['r_proj'], bins=100, color='salmon', edgecolor='black', alpha=0.7)
    axes[0, 1].set_title(r'Radial Projection Distribution ($\Delta R$)')
    axes[0, 1].set_xlabel(r'$R_{proj}$')
    axes[0, 1].set_yscale('log')
    axes[0, 1].grid(True, alpha=0.3)

    # 3. 切向偏移分布
    axes[1, 0].hist(res['t_proj'], bins=100, color='lightgreen', edgecolor='black', alpha=0.7)
    axes[1, 0].set_title('Transverse Projection Distribution ($\Theta_R$)')
    axes[1, 0].set_xlabel(r'$T_{proj}$')
    axes[1, 0].set_yscale('log')
    axes[1, 0].grid(True, alpha=0.3)

    # 4. 夹角分布
    axes[1, 1].hist(res['theta'], bins=90, range=(0, 180), color='plum', edgecolor='black', alpha=0.7)
    axes[1, 1].set_title(r'Angle Distribution ($\phi$)')
    axes[1, 1].set_xlabel(r'Angle [Degrees]')
    axes[1, 1].set_xticks(np.arange(0, 181, 30))
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
else:
    print("No comparison data to plot.")

# %%
plot_field_from_file('/mnt/18T/output_2D/lc_c10/3000_1024_onerun/0.000_delta_c.bin', title='onerun')
plot_field_from_file('/mnt/18T/output_2D/lc_c10/3000_1024_runtime/lightcone_post_m1_delta.bin', title='runtime')
# %%
header, res = load_lc_comparison('/mnt/18T/output_2D/lc_c10/comparsion/runtime_runtime_2_one_run.bin')
if res is not None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=200)
    fig.suptitle(f"Lightcone Comparison Statistics\n(Match: {header['match_count']}, Missing A: {header['mis_b']}, Missing B: {header['mis_a']}) \n Runtime .vs. OneRun ", fontsize=14)

    # 1. 绝对距离分布
    axes[0, 0].hist(res['dist'], bins=100, color='skyblue', edgecolor='black', alpha=0.7)
    axes[0, 0].axvline(np.mean(res['dist']), color='red', linestyle='--', alpha=0.5)
    axes[0, 0].set_title(r'Absolute Distance Distribution (Grid Units)')
    axes[0, 0].set_xlabel(r'$\Delta r$')
    axes[0, 0].set_yscale('log')
    axes[0, 0].grid(True, alpha=0.3)

    # 2. 径向偏移分布
    axes[0, 1].hist(res['r_proj'], bins=100, color='salmon', edgecolor='black', alpha=0.7)
    axes[0, 1].axvline(np.mean(res['r_proj']), color='red', linestyle='--', alpha=0.5)
    axes[0, 1].set_title(r'Radial Projection Distribution ($\Delta R$)')
    axes[0, 1].set_xlabel(r'$R_{proj}$')
    axes[0, 1].set_yscale('log')
    axes[0, 1].grid(True, alpha=0.3)

    # 3. 切向偏移分布
    axes[1, 0].hist(res['t_proj'], bins=100, color='lightgreen', edgecolor='black', alpha=0.7)
    axes[1, 0].axvline(np.mean(res['t_proj']), color='red', linestyle='--', alpha=0.5)
    axes[1, 0].set_title('Transverse Projection Distribution ($\Theta_R$)')
    axes[1, 0].set_xlabel(r'$T_{proj}$')
    axes[1, 0].set_yscale('log')
    axes[1, 0].grid(True, alpha=0.3)

    # 4. 夹角分布
    axes[1, 1].hist(res['theta'], bins=90, range=(0, 180), color='plum', edgecolor='black', alpha=0.7)
    axes[1, 1].axvline(np.mean(res['theta']), color='red', linestyle='--', alpha=0.5)
    axes[1, 1].set_title(r'Angle Distribution ($\phi$)')
    axes[1, 1].set_xlabel(r'Angle [Degrees]')
    axes[1, 1].set_xticks(np.arange(0, 181, 30))
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
else:
    print("No comparison data to plot.")


# %% 

# 绘制生长率比值
end = 50
k0, xi0 = loadpower(f'{Path}{z_list[-1]:.3f}_power.bin')
growth_ratio = np.zeros_like(z_list)
for i, z in enumerate(z_list):
    k, xi = loadpower(f'{Path}{z:.3f}_power.bin')
    growth_ratio[i] = np.mean(xi[1:end, 2] / xi0[1:end, 2] )

def growth_ratio_cube(a):
    om = 0.2606667599+0.0489746816
    ol = 1- om
    hsq=om/a**3+(1-om-ol)/a**2+ol
    oma=om/(a**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
    Dgrow=a*ga/g
    return Dgrow

growth_ratio_3D = growth_ratio_cube(1/(1+z_list))
plt.plot(z_list, growth_ratio)
plt.plot(z_list, growth_ratio_3D, linestyle='--')
plt.legend(['3D', 'Cube'])
plt.show()

# %%
