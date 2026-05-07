# 2D纯CDM宇宙学模型功率谱生成
# 在2D纯CDM宇宙中：
# - 膨胀律：H(a) = H0 * a^-1（因为密度∝a^-2）
# - 增长因子：D(a) ∝ a（线性增长）
# - Dgrow_ratio = D(a2)/D(a1) = a2/a1
# 初始功率谱 = P(k, z=0) * Dgrow_ratio^2

import camb,sys,os,re,time
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d

test = 0
all_reps = False
all_camb = True
n=50
ni=2

def match_para(para):
    file_path = './parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()
    pattern = r'parameter\s*::\s*%s\s*=\s*([^\s]+)'%para
    match = re.search(pattern, content)
    if match:
        variable_value = match.group(1).strip()
        return float(variable_value)
    else:
        print(para,"Pattern not found in the file.")
        sys.exit()

def get_z(n):
    z = [0]*n
    a = [0]*n
    z[0] = z_max
    a = np.linspace(1/(1+z[0]),1,n+1)
    z=1/a-1
    return z[::-1]

def get_Pk_nonlin_CDM_2D(n):
    Pk_nonlin_CDM = [0]*len(z_nonlin)
    nz = int(np.ceil(len(z_nonlin)/n_PK))
    for i in range(0,nz):
        z_n=z_nonlin[i*n_PK:(i+1)*n_PK]
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, 
                          neutrino_hierarchy='degenerate', num_massive_neutrinos=0, 
                          mnu=0, nnu=3.044, standard_neutrino_neff=0)
        pars.InitPower.set_params(As=As,ns=ns)
        print('i = %d/%d;z_max = %.3f;len = %d'%(i+1,nz,z_n[-1],len(z_n)))
        pars.set_matter_power(redshifts=z_n, kmax=k_ic_max, nonlinear=True)
        pars.NonLinear = camb.model.NonLinear_both
        results = camb.get_results(pars)
        kh_nonlin, z_nonlin[i*n_PK:i*n_PK+len(z_n)],Pk_nonlin_CDM[i*n_PK:i*n_PK+len(z_n)]= \
            results.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,
                                              var1='delta_nonu',var2='delta_nonu')
        pars = 0
    return Pk_nonlin_CDM[::-1],kh_nonlin,z_nonlin[::-1]

# 读取参数
if (1):
    H0=match_para('h0')*100
    omega_bar=match_para('omega_bar')
    omega_cdm=match_para('omega_cdm')
    omega_r=match_para('omega_r')
    omk=0.0

    m_nus = np.array([0.0, 0.0, 0.0])
    mnu = 0.0
    neutrino_hierarchy = 'degenerate'
    num_massive_neutrinos = 0
    nnu = 3.044
    standard_neutrino_neff = 0
    ns = match_para('n_s')
    As = match_para('A_s')
    z_max = np.atleast_1d(np.loadtxt('./z_checkpoint.txt'))[0]

    pi = np.pi
    ratio_cs = match_para('ratio_cs')
    ng = match_para('ng')
    box = match_para('box')
    nnt = match_para('nnt')
    # nn = match_para('nn')
    ngp = ng/nnt
    ngb = match_para('ngb')
    istep_max = match_para('istep_max')
    ngt = ngp+2*ngb 
    npbin  = 400
    kmin = 1e-4
    kmax = 20
    ombh2 = omega_bar*(H0/100)**2
    omch2 = omega_cdm*(H0/100)**2

    sigma_nu = (4/11)**(1/3)
    N_nu = 3
    N_eff = 0
    k_b = 8.617342e-5
    T_gama = 2.7255
    Mass_nu = 0.0
    omega_m = omega_bar+omega_cdm
    omega_cb = omega_bar+omega_cdm
    omega_l = 1-omega_bar-omega_cdm-omega_r
    f_nu = 0.0
    
    # 2D纯CDM宇宙学参数
    # H(a) = H0 * a^-1
    # D(a) = a (线性增长)
    # omHsq_2D = dt/da * a^2 * H(a) = dt/da * a * H0
    # 对于纯CDM: da/dt = H0 * a^-1, 所以 dt/da = a/H0
    # omHsq_2D = (a/H0) * a * H0 = a^2
    # 但我们需要一个常数来控制时间步长
    omHsq_2D = 1.0/H0

    file_path = './parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()
    pattern = r'parameter\s*::\s*opath\s*=\s*([^\s]+)'
    match = re.search(pattern, content)
    opath = os.path.expanduser(match.group(1).strip()[1:-1])
    print(f"Variable value of opath : {opath}")
    nupath = opath+"neutrinos"
    print(f"Variable value of nupath : {nupath}")
    print('\n'+('+'*40+'\n')*2)
    print('2D Pure CDM Cosmology Paras:\n\n   omega_r:   %.6f\n   omega_b:   %.6f\n   omega_c:   %.6f\n   omega_l:   %.6f\n   mass_nu:   %.3f              eV\n      f_nu:   %.6f\n\n'%(omega_r,omega_bar,omega_cdm,omega_l,Mass_nu,f_nu))
    print('Simulation Paras:\n\n     npbin:   %d\n      kmin:   %.3f\n      kmax:   %.3f\n\n\n'%(npbin,kmin,kmax))

    kh_nonlin = np.exp(np.linspace(np.log(6*kmin), np.log(kmax),  npbin))

# 创建目录
try:
    os.system('mkdir -p '+nupath+'/TF')
    os.system('mkdir -p '+nupath+'/Pk_m')
    os.system('mkdir -p '+nupath+'/Pk_nu')
    os.system('mkdir -p '+nupath+'/Pk_nus')
    os.system('mkdir -p %s/IC'%(nupath))
except:
    None

# 设置红移数组
if test:
    n=n*ni
else:
    n=n
n_PK = 150
k_ic_min = 1e-4
k_ic_max = max(1e2,kmax*1.5)

z_nonlin = get_z(n)

# 获取红移0的功率谱（用CAMB）
print('Getting P(k) at z=0 from CAMB...')
par = camb.CAMBparams()
par.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, 
                  neutrino_hierarchy='degenerate', num_massive_neutrinos=0, 
                  mnu=0, nnu=3.044, standard_neutrino_neff=0)
par.InitPower.set_params(As=As,ns=ns)
par.set_matter_power(redshifts=[0], kmax=k_ic_max, nonlinear=True)
par.NonLinear = camb.model.NonLinear_both
result = camb.get_results(par)
kh_ic, _, Pk_cb_z0 = result.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, 
                                                       npoints=npbin, var1='delta_nonu', var2='delta_nonu')
Pk_cb_z0 = Pk_cb_z0[0]
print(f'P(k) at z=0 obtained, k range: {kh_ic.min():.4f} - {kh_ic.max():.4f} h/Mpc')

# 2D纯CDM宇宙学函数
# H(a) = H0 * a^-1
def Ha_2D(a):
    return H0 / a

# 增长因子 D(a) = a（线性增长）
def Dgrow_2D(a):
    return a

# 增长因子比值 D(a2)/D(a1) = a2/a1
def Dgrow_ratio_2D(a1, a2):
    return a2 / a1

# 共形时间 tau(a) = ∫ da / (a^2 * H(a)) = ∫ da / (a * H0) = ln(a) / H0
def taua_2D(a):
    return np.log(a) / H0

# 共动距离 chi(a) = c * ∫ da / (a^2 * H(a)) = c * ln(a) / H0
def chia_2D(a):
    c = 299792.458  # km/s
    return c * np.log(a) / H0

# 计算膨胀历史（使用2D纯CDM模型）
print('Calculating Expansion History for 2D pure CDM...')
n_a = int(istep_max)
omHsq = 1/H0# 2/3*np.sqrt(1/omega_m)/H0
dt0 = 2.5e-3
t = -np.arange(n_a)*dt0
c2h = 2997.92458  # 光速，单位 100km/s
chi_ex = 1/np.zeros(n_a)
tau = np.zeros(n_a)
a_ex = np.zeros(n_a)
a_ex[0] = 1
tau[0] = taua_2D(1)
chi_ex[0] = 0.0

T1 = time.time()
for i in range(n_a-1):
    Hai = Ha_2D(a_ex[i])
    a_ex[i+1] = a_ex[i] - omHsq*Hai*a_ex[i]**3 * dt0
    
    # 共形时间
    tau[i+1] = taua_2D(a_ex[i+1])
    
    # 共动距离 chi，使用公式 dchi = c/a * dt
    chi_ex[i+1] = chi_ex[i] + c2h / a_ex[i] * dt0
    
    if (a_ex[i+1] <= 1./301 or np.isnan(a_ex[i+1])):
        break

i_end = i
T2 = time.time()
print("EH:\n\n      time:   %.2f seconds\n      step:   %d\n"%(T2-T1,i_end))
print(a_ex.max(),a_ex.min())
print(tau.max(),tau.min())
print(chi_ex.max(),chi_ex.min())
idx =  np.where(a_ex<1/2)[0][0]
print(idx)
print(a_ex[idx])
print(chi_ex[idx])



# 保存膨胀历史
os.system('rm '+nupath+'/*.txt 2>/dev/null')
np.savetxt(nupath+'/s_a_tau_H.txt',np.array([t,a_ex,tau,chi_ex]))
print("      save:   '%s'\n\n\n"%nupath+'/s_a_tau_H.txt')


# 保存功率谱
z_powerpoint = open(nupath+'/z_powerpoint.txt', 'w')
z_values = open(nupath+'/z_values.txt', 'w')
k_values = open(nupath+'/k_values.txt', 'w')

# 保存z=0的功率谱作为初始功率谱
Pk_cb_ic = Pk_cb_z0 * Dgrow_2D(1/(1+z_max))**2
np.savetxt(nupath+'/IC/Pcb_ic.txt', np.array([kh_ic, Pk_cb_ic]).T)
# 对于纯CDM，中微子功率谱为0
np.savetxt(nupath+'/IC/Pnu_ic.txt', np.zeros(len(kh_nonlin)))

for i in range(len(kh_nonlin)):
    k_values.write('%3.12f\n'%kh_nonlin[i])

if mnu>0:
    # 计算各红移的功率谱（使用2D纯CDM增长因子）
    print('Calculating P(k) at different redshifts using 2D pure CDM growth factor...')
    n = int(z_nonlin.shape[0])
    Pk_nl = [0]*n

    # 插值z=0的功率谱
    Pk_z0_interp = interp1d(kh_ic, Pk_cb_z0, kind='linear', bounds_error=False, fill_value="extrapolate")

    for i in range(n):
        z_i = z_nonlin[i]
        a_i = 1.0 / (1.0 + z_i)
        
        # 2D纯CDM: P(k, z) = P(k, z=0) * D(a)^2 = P(k, z=0) * a^2
        D2 = Dgrow_2D(a_i)**2
        Pk_nl[i] = Pk_z0_interp(kh_nonlin) * D2
        
        if np.any(Pk_nl[i] < 0):
            print('Pk has negative values at z =', z_i)
            exit()

    print(f'P(k) calculated for {n} redshifts')
    z_str = ''
    for i in range(n):
        np.savetxt(nupath+'/Pk_cb_%3.4f.txt'%z_nonlin[i], Pk_nl[i])
        
        if np.any(Pk_nl[i] < 0):
            print('Pk_cb has negative values')
            print(Pk_nl[i])
            print(z_nonlin[i])
            exit()
        
        # 纯CDM没有中微子
        np.savetxt(nupath+'/Pk_nu_%3.4f.txt'%z_nonlin[i], np.zeros(len(kh_nonlin)))
        # 传递函数为1
        np.savetxt(nupath+'/Tf_nu_%3.4f.txt'%z_nonlin[i], np.ones(len(kh_nonlin)))
        z_str += '%3.4f\n'%z_nonlin[i]
    z_powerpoint.write(z_str)

z_powerpoint.close()
z_values.close()
k_values.close()

print('******************** Pk_2D done ********************\n\n\n')
