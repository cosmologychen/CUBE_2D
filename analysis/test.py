#%% 无锡答辩
# 无锡答辩
import matplotlib.pyplot as plt
import pandas as pd

# 1. 按照要求重构数据字典: {Name: [Year, Np, Mp, Type]}
# 注：Mp 均为 CDM 粒子质量
data = {
    'Our Plan': [2026, 12288**3, 1.251e6, 'Our Plan'],
    'Uchuu': [2020, 12800**3,3.7e8, 'Other Works'],
    # 'Shin-Uchuu': [2020, 600**3, 8.97e5, 'Other Works'],
    # 'FLAMINGO': [2023, 5040**3, 1.0e9, 'Hydro'],
    'TianNu': [2017, 2.97e12, 4.8e7, 'Other Works'],
    'Far Point':[2021, 1.86e12, 4.6e7, 'Other Works'],
    'Cosmo-pi': [2019, 4.39e12, 6.2e8, 'Other Works'],
    # 'MillenniumTNG': [2023, 320**3, 1.32e8, 'Hydro'],
    'Euclid Flagship 1': [2017, 12600**3, 2.4e9, 'Other Works'],
    'Euclid Flagship 2': [2025, 4.0e12, 1e9, 'Other Works'],
    'Last Journey': [2020, 1.24e12, 2.7e9, 'Other Works'],
    'OuterRim ': [2019, 1.07e12, 1.85e9, 'Other Works'],
    'Dark Sky': [2014, 1.1e12, 3.8e10, 'Other Works']
}

# 转换为 DataFrame 方便绘图
df = pd.DataFrame.from_dict(data, orient='index', columns=['Year', 'Np', 'Mp', 'Type']).reset_index()
df.rename(columns={'index': 'Name'}, inplace=True)

# 2. 绘图设置
plt.figure(figsize=(13, 9), dpi=150)
colors = {'Other Works': '#1f77b4', 'Hydro': '#2ca02c','Our Plan':'#ff7f0e'}
markers = {'Other Works': 'o', 'Hydro': 'o', 'Our Plan': '*'}
sizes = {'Other Works': 250, 'Hydro': 250, 'Our Plan': 1500}

# 3. 绘制散点
for t in df['Type'].unique():
    mask = df['Type'] == t
    plt.scatter(df[mask]['Np'], df[mask]['Mp'], s=sizes[t], label=t, marker=markers[t],
                color=colors[t], alpha=0.8, edgecolors='white', linewidth=1.5, zorder=3)

# 4. 手动处理标签避让（解决 Far Point & TianNu 重叠）
offsets = {
    'Far Point': (-12, -12), 
    'TianNu': (10, 5),
    'Uchuu': (1, -30),
    'Sahyadri': (0, -20),
    'COLIBRE (m6)': (0, 15),
    'Cosmo-pi': (-10, -20),
    'Euclid Flagship 2': (-10, -10),
    'OuterRim ': (0, -30),
    'Our Plan': (0, 20),
    'Last Journey': (0, 15),
    'Euclid Flagship 1': (15, 0),
}

for i, row in df.iterrows():
    name = row['Name']
    x_off, y_off = offsets.get(name, (10, 5))
    plt.annotate(
        f"{name} ({int(row['Year'])})", 
        (row['Np'], row['Mp']),
        textcoords="offset points", 
        xytext=(x_off, y_off),
        ha='left' if x_off >= 0 else 'right',
        fontsize=24, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0, ec='none'),
        zorder=4
    )

# 5. 图形修饰
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

plt.xlabel('Total Particle Number ($N_p$)', fontsize=24)
plt.ylabel('CDM Particle Mass ($M_{CDM}$) [$M_{\odot}/h$]', fontsize=24)

plt.grid(True, which="both", ls="--", alpha=0.3, zorder=0)
plt.legend( loc='upper right', frameon=True, shadow=True, fontsize=24)

plt.tight_layout()
plt.show()
# %%
