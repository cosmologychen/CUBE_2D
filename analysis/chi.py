# %%
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# 设置路径和观测者红移列表
Path = '/mnt/18T/output_2D/lc_c10/3000_1024_iclc/'
z_list = [0, 1, 5, 10, 100]

# 读取膨胀历史数据
eh_file = os.path.join(Path, 'neutrinos/s_a_tau_H.txt')
eh_data = np.loadtxt(eh_file)
eh_t, eh_a, eh_eta, eh_chi, eh_D = eh_data

# 1. 数据预处理
mask = (eh_a > 0) & (~np.isnan(eh_a)) & (eh_a <= 1.0)
eh_t, eh_a, eh_eta, eh_chi, eh_D = eh_t[mask], eh_a[mask], eh_eta[mask], eh_chi[mask], eh_D[mask]

eh_z = 1.0 / eh_a - 1.0

# 2. 建立插值函数
z_to_chi = interp1d(eh_z, eh_chi, kind='linear', fill_value='extrapolate')
chi_to_t = interp1d(eh_chi, eh_t, kind='linear', fill_value='extrapolate')
chi_to_a = interp1d(eh_chi, eh_a, kind='linear', fill_value='extrapolate')
chi_to_eta = interp1d(eh_chi, eh_eta, kind='linear', fill_value='extrapolate')
chi_to_D = interp1d(eh_chi, eh_D, kind='linear', fill_value='extrapolate')

# 3. 绘图计算
r = np.linspace(0, 3000, 500)

fig, axes = plt.subplots(3, 1, figsize=(6, 20))
titles = [r'Drowth Factor $D/D_{obs}$', r'Scale Factor $a/a_{obs}$', 'Conformal Time $\eta$']
ylabels = [r'$D/D_{obs}$', r'$a/a_{obs}$', '$\eta$']

for z_obs in z_list:
    try:
        chi_obs = z_to_chi(z_obs)
    except:
        continue
        
    chi_look = chi_obs + r
    
    # 只取在数据表范围内的部分
    valid = (chi_look >= eh_chi.min()) & (chi_look <= eh_chi.max())
    if not np.any(valid):
        continue
        
    curr_r = r[valid]
    curr_chi = chi_look[valid]
    
    # 提取光锥上的物理量
    D_vals = chi_to_D(curr_chi)
    a_vals = chi_to_a(curr_chi)
    eta_vals = chi_to_eta(curr_chi)
    
    label = f'$z_{{obs}}={z_obs}$'
    axes[0].plot(curr_r, D_vals/D_vals[0], label=label)
    axes[1].plot(curr_r, a_vals*(1+z_obs), label=label)
    axes[2].plot(curr_r, eta_vals, label=label)
    axes[0].set_ylim(0,1.1)
    axes[1].set_yscale('log')
    

# 4. 修饰图表
for i in range(3):
    axes[i].set_title(titles[i])
    axes[i].set_xlabel('Relative Comoving Distance $r$ (Mpc/h)')
    axes[i].set_ylabel(ylabels[i])
    axes[i].legend()
    axes[i].grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()

# %%

# 读取膨胀历史数据
eh_file = os.path.join(Path, 'neutrinos/s_a_tau_H.txt')
eh_data = np.loadtxt(eh_file)
eh_t, eh_a, eh_eta, eh_chi, eh_D = eh_data

# 1. 数据预处理
mask = (eh_a > 0) & (~np.isnan(eh_a)) & (eh_a <= 1.0)
eh_t, eh_a, eh_eta, eh_chi, eh_D = eh_t[mask], eh_a[mask], eh_eta[mask], eh_chi[mask], eh_D[mask]
chi_to_t = interp1d(eh_chi, eh_t, kind='linear', fill_value='extrapolate')
a_to_t = interp1d(eh_a, eh_t, kind='linear', fill_value='extrapolate')
t_lc = np.zeros_like(r)
for i in range(len(r)):
    t_lc[i] = chi_to_t(r[i])
    
for z_obs in z_list:
    t = a_to_t(1/(1+z_obs))

    
    label = f'$z_{{obs}}={z_obs}$'
    plt.plot(r, t_lc+t, label=label)

plt.title(r'Super-Conformal Time $\tau$')
plt.xlabel('Relative Comoving Distance $r$ (Mpc/h)')
plt.ylabel(r'$\tau$')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
# %%
chi_to_t = interp1d(eh_chi, eh_t, kind='linear', fill_value='extrapolate')
a_to_t = interp1d(eh_a, eh_t, kind='linear', fill_value='extrapolate')
t_to_a = interp1d(eh_t, eh_a, kind='linear', fill_value='extrapolate')
t_lc = np.zeros_like(r)
for i in range(len(r)):
    t_lc[i] = chi_to_t(r[i])
    
for z_obs in z_list:
    t = a_to_t(1/(1+z_obs))
    
    label = f'$z_{{obs}}={z_obs}$'
    plt.plot(r, t_to_a((t_lc+t))*(1+z_obs), label=label)

plt.title(r'Scale Factor $a/a_{obs}$')
plt.xlabel('Relative Comoving Distance $r$ (Mpc/h)')
plt.ylabel(r'$a/a_{obs}$')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
 
# %%
