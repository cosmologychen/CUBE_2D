# %% lightcone_analysis.py - Lightcone粒子分析
# lightcone_analysis.py - Lightcone粒子分析
# 
# 功能：
# 1. 画图对比两个lightcone粒子列表的位置和密度场
# 2. 计算两个lightcone粒子列表对应粒子的位移矢量，分解为径向和切向
# 3. 统计dis, dr, theta三个距离的分布和互相关

import numpy as np
import matplotlib.pyplot as plt
import os
import struct
from scipy.stats import histogram

# %% 设置绘图参数
# 设置绘图参数
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
    'axes.spines.right': False,
    'axes.spines.top': False,
    'axes.grid': False,
    'axes.labelpad': 8.0,
    'axes.labelsize': 10,
    'axes.labelcolor': 'k',
    'axes.axisbelow': True,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'font.size': 10,
    'legend.fontsize': 12,
    'legend.frameon': False,
    'lines.markersize': 6,
})

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

# %% 工具函数
# 工具函数

def loadfield2d(fn):
    """
    加载2D场
    fn: filename
    Returns: a, n (场数据, 尺寸)
    """
    fid = open(fn, 'rb')
    p1 = np.fromfile(fid, dtype=np.float32)
    fid.close()
    n = round(len(p1) ** (1/2))
    if (n**2 != len(p1)):
        print('shape no match')
        print(n, n**2, len(p1))
        return 0, 0
    a = np.reshape(p1, (n, n), order='F')
    return a, n

def load2dxpos(fn):
    """
    加载2D粒子位置
    fn: filename
    Returns: pos, n (位置数组, 粒子数开方)
    """
    fid = open(fn, 'rb')
    p1 = np.fromfile(fid, dtype=np.float32)
    fid.close()
    n = round(len(p1)/2)
    if (n*2 != len(p1)):
        print('shape no match')
        print(n, n**2, len(p1))
        return 0, 0
    a = np.reshape(p1, (2, n), order='F')
    return a, int(np.sqrt(n))

def load_pid(fn):
    """
    加载粒子ID
    fn: filename
    Returns: pid (粒子ID数组)
    """
    fid = open(fn, 'rb')
    pid = np.fromfile(fid, dtype=np.int64)
    fid.close()
    return pid

def get_sim_info(prefix):
    """
    获取模拟信息
    prefix: 文件前缀
    Returns: sim (模拟信息字典)
    """
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
    sim['ng_global'] = sim['nnt'] * sim['nt'] * 4
    return sim

# %% 主分析类
# 主分析类

class LightconeAnalysis:
    def __init__(self, path1, path2=None, observer=(0, 0)):
        """
        初始化分析
        path1: 第一个lightcone数据路径
        path2: 第二个lightcone数据路径（可选）
        observer: 观测者位置
        """
        self.path1 = path1
        self.path2 = path2
        self.observer = observer
        
        # 加载数据
        self.load_data()
    
    def load_data(self):
        """加载lightcone数据"""
        print('Loading lightcone data...')
        
        # 加载第一个lightcone
        self.pos1, n1 = load2dxpos(self.path1 + '_xp.bin')
        self.vp1, _ = load2dxpos(self.path1 + '_vp.bin')
        try:
            self.pid1 = load_pid(self.path1 + '_pid.bin')
        except:
            self.pid1 = np.arange(self.pos1.shape[1])
        print(f'  Lightcone 1: {self.pos1.shape[1]} particles')
        
        # 加载第二个lightcone（如果存在）
        if self.path2:
            self.pos2, n2 = load2dxpos(self.path2 + '_xp.bin')
            self.vp2, _ = load2dxpos(self.path2 + '_vp.bin')
            try:
                self.pid2 = load_pid(self.path2 + '_pid.bin')
            except:
                self.pid2 = np.arange(self.pos2.shape[1])
            print(f'  Lightcone 2: {self.pos2.shape[1]} particles')
    
    def plot_positions(self, save_path=None):
        """
        画图对比两个lightcone粒子列表的位置
        save_path: 保存路径
        """
        print('Plotting positions...')
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 6), dpi=150)
        
        # 第一个lightcone
        ax1 = axes[0]
        ax1.scatter(self.pos1[0, :], self.pos1[1, :], c='black', s=0.1, marker='o', 
                   edgecolors='none', alpha=0.5)
        ax1.scatter([self.observer[0]], [self.observer[1]], c='red', s=50, marker='x', 
                   label='Observer')
        ax1.set_aspect('equal')
        ax1.set_title('Lightcone 1')
        ax1.set_xlabel('x [Mpc/h]')
        ax1.set_ylabel('y [Mpc/h]')
        ax1.legend()
        
        # 第二个lightcone
        if self.path2:
            ax2 = axes[1]
            ax2.scatter(self.pos2[0, :], self.pos2[1, :], c='blue', s=0.1, marker='o', 
                       edgecolors='none', alpha=0.5)
            ax2.scatter([self.observer[0]], [self.observer[1]], c='red', s=50, marker='x', 
                       label='Observer')
            ax2.set_aspect('equal')
            ax2.set_title('Lightcone 2')
            ax2.set_xlabel('x [Mpc/h]')
            ax2.set_ylabel('y [Mpc/h]')
            ax2.legend()
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
            print(f'  Saved to {save_path}')
        else:
            plt.show()
        
        plt.close()
    
    def plot_density_field(self, grid_size=256, box_size=None, save_path=None):
        """
        画图对比两个lightcone粒子列表的密度场
        grid_size: 网格大小
        box_size: 盒子大小
        save_path: 保存路径
        """
        print('Plotting density fields...')
        
        if box_size is None:
            box_size = max(np.max(self.pos1), np.max(self.pos2 if self.path2 else self.pos1)) * 1.1
        
        fig, axes = plt.subplots(1, 2 if self.path2 else 1, figsize=(12 if self.path2 else 6, 6), dpi=150)
        if not self.path2:
            axes = [axes]
        
        # 计算密度场
        for i, (pos, ax, title) in enumerate([
            (self.pos1, axes[0], 'Lightcone 1'),
            *([(self.pos2, axes[1], 'Lightcone 2')] if self.path2 else [])
        ]):
            # 创建网格
            density = np.zeros((grid_size, grid_size))
            dx = box_size / grid_size
            
            # 统计每个格点的粒子数
            x_idx = np.clip((pos[0, :] / box_size * grid_size).astype(int), 0, grid_size-1)
            y_idx = np.clip((pos[1, :] / box_size * grid_size).astype(int), 0, grid_size-1)
            
            for j in range(len(x_idx)):
                density[x_idx[j], y_idx[j]] += 1
            
            # 计算密度扰动
            mean_density = len(x_idx) / (grid_size ** 2)
            delta = (density - mean_density) / mean_density
            
            # 绘图
            im = ax.imshow(delta.T, origin='lower', extent=[0, box_size, 0, box_size],
                          cmap='viridis', vmin=-1, vmax=3)
            ax.set_aspect('equal')
            ax.set_title(f'{title} Density Field')
            ax.set_xlabel('x [Mpc/h]')
            ax.set_ylabel('y [Mpc/h]')
            plt.colorbar(im, ax=ax, label=r'$\delta$')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
            print(f'  Saved to {save_path}')
        else:
            plt.show()
        
        plt.close()
    
    def calculate_displacement(self):
        """
        计算两个lightcone粒子列表对应粒子的位移矢量
        并分解为径向和切向两个方向
        
        Returns: dis, dr, theta, R (总位移, 径向位移, 切向位移, 粒子距离)
        """
        print('Calculating displacement vectors...')
        
        if not self.path2:
            print('  Error: Need two lightcone datasets')
            return None, None, None, None
        
        # 找到共同的粒子ID
        common_pid = np.intersect1d(self.pid1, self.pid2)
        print(f'  Found {len(common_pid)} common particles')
        
        if len(common_pid) == 0:
            print('  Error: No common particles found')
            return None, None, None, None
        
        # 获取共同粒子的索引
        idx1 = np.searchsorted(self.pid1, common_pid)
        idx2 = np.searchsorted(self.pid2, common_pid)
        
        # 计算位移矢量
        dx = self.pos2[0, idx2] - self.pos1[0, idx1]
        dy = self.pos2[1, idx2] - self.pos1[1, idx1]
        
        # 计算粒子到观测者的距离（使用第一个lightcone的位置）
        R = np.sqrt((self.pos1[0, idx1] - self.observer[0])**2 + 
                    (self.pos1[1, idx1] - self.observer[1])**2)
        
        # 计算径向单位向量
        r_hat_x = (self.pos1[0, idx1] - self.observer[0]) / (R + 1e-10)
        r_hat_y = (self.pos1[1, idx1] - self.observer[1]) / (R + 1e-10)
        
        # 切向单位向量（垂直于径向）
        t_hat_x = -r_hat_y
        t_hat_y = r_hat_x
        
        # 分解位移
        dr = dx * r_hat_x + dy * r_hat_y  # 径向位移
        dt = dx * t_hat_x + dy * t_hat_y  # 切向位移
        
        # 总位移
        dis = np.sqrt(dx**2 + dy**2)
        
        # 切向角度 theta = dt / R
        theta = np.abs(dt) / (R + 1e-10)
        
        print(f'  Mean displacement: {np.mean(dis):.4f} Mpc/h')
        print(f'  Mean radial displacement: {np.mean(dr):.4f} Mpc/h')
        print(f'  Mean tangential angle: {np.mean(theta):.6f} rad')
        
        return dis, dr, theta, R
    
    def plot_statistics(self, dis, dr, theta, R, save_path=None):
        """
        统计dis, dr, theta三个距离的分布和互相关
        """
        print('Plotting statistics...')
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10), dpi=150)
        
        # dis分布
        ax = axes[0, 0]
        ax.hist(dis, bins=100, density=True, alpha=0.7, color=colors[0])
        ax.set_xlabel('dis [Mpc/h]')
        ax.set_ylabel('PDF')
        ax.set_title('Total Displacement Distribution')
        ax.set_yscale('log')
        
        # dr分布
        ax = axes[0, 1]
        ax.hist(dr, bins=100, density=True, alpha=0.7, color=colors[1])
        ax.set_xlabel('dr [Mpc/h]')
        ax.set_ylabel('PDF')
        ax.set_title('Radial Displacement Distribution')
        
        # theta分布
        ax = axes[0, 2]
        ax.hist(theta, bins=100, density=True, alpha=0.7, color=colors[2])
        ax.set_xlabel('theta [rad]')
        ax.set_ylabel('PDF')
        ax.set_title('Tangential Angle Distribution')
        
        # dis vs R
        ax = axes[1, 0]
        ax.scatter(R, dis, s=0.1, alpha=0.5, c=colors[0])
        ax.set_xlabel('R [Mpc/h]')
        ax.set_ylabel('dis [Mpc/h]')
        ax.set_title('dis vs R')
        
        # dr vs R
        ax = axes[1, 1]
        ax.scatter(R, dr, s=0.1, alpha=0.5, c=colors[1])
        ax.set_xlabel('R [Mpc/h]')
        ax.set_ylabel('dr [Mpc/h]')
        ax.set_title('dr vs R')
        
        # theta vs R
        ax = axes[1, 2]
        ax.scatter(R, theta, s=0.1, alpha=0.5, c=colors[2])
        ax.set_xlabel('R [Mpc/h]')
        ax.set_ylabel('theta [rad]')
        ax.set_title('theta vs R')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
            print(f'  Saved to {save_path}')
        else:
            plt.show()
        
        plt.close()
        
        # 计算互相关
        print('\nCorrelation Analysis:')
        print(f'  corr(dis, R) = {np.corrcoef(dis, R)[0, 1]:.4f}')
        print(f'  corr(dr, R) = {np.corrcoef(dr, R)[0, 1]:.4f}')
        print(f'  corr(theta, R) = {np.corrcoef(theta, R)[0, 1]:.4f}')
        print(f'  corr(dis, dr) = {np.corrcoef(dis, dr)[0, 1]:.4f}')
        print(f'  corr(dis, theta) = {np.corrcoef(dis, theta)[0, 1]:.4f}')
        print(f'  corr(dr, theta) = {np.corrcoef(dr, theta)[0, 1]:.4f}')

# %% 主程序
# 主程序

if __name__ == '__main__':
    # 示例用法
    # 修改以下路径为实际路径
    
    # Path = '/mnt/18T/output_2D/600_3072/'
    # 
    # # 创建分析对象
    # analysis = LightconeAnalysis(
    #     path1=Path + 'lightcone_1',
    #     path2=Path + 'lightcone_2',
    #     observer=(Path_info['box']/2, Path_info['box']/2)
    # )
    # 
    # # 绘制位置对比
    # analysis.plot_positions(save_path=Path + 'fig/lightcone_positions.pdf')
    # 
    # # 绘制密度场对比
    # analysis.plot_density_field(save_path=Path + 'fig/lightcone_density.pdf')
    # 
    # # 计算位移
    # dis, dr, theta, R = analysis.calculate_displacement()
    # 
    # # 绘制统计
    # analysis.plot_statistics(dis, dr, theta, R, save_path=Path + 'fig/lightcone_statistics.pdf')
    
    print('Lightcone Analysis Module Loaded')
    print('Usage:')
    print('  analysis = LightconeAnalysis(path1, path2, observer=(x, y))')
    print('  analysis.plot_positions()')
    print('  analysis.plot_density_field()')
    print('  dis, dr, theta, R = analysis.calculate_displacement()')
    print('  analysis.plot_statistics(dis, dr, theta, R)')
