import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

from DataProcess import SignalAnalysis, density_string_integral

# Load data
shot = '232094'
data = sio.loadmat(shot+'.mat')

Ip =data['ip1'].reshape(1,-1)[0]
t_Ip = data['t_ip1'].reshape(1,-1)[0]


# Get data
sig1 = data['sig1'].reshape(1,-1)[0]
sig2 = data['sig2'].reshape(1,-1)[0]
sig3 = data['sig3'].reshape(1,-1)[0]
ref = data['ref'].reshape(1,-1)[0]


lendata = len(sig1)
time_tot = 100*1e-3 # Total time,s
time_step = time_tot/lendata # Time step,s

fs        = 60   *1e6
fm_target = 2.125*1e6


SA = SignalAnalysis(
    fs              = fs,
    fm_target       = fm_target,
    signal_plasma   = sig1,
    signal_ref      = ref
)

SA2 = SignalAnalysis(
    fs              = fs,
    fm_target       = fm_target,
    signal_plasma   = sig2,
    signal_ref      = ref
)

SA3 = SignalAnalysis(
    fs              = fs,
    fm_target       = fm_target,
    signal_plasma   = sig3,
    signal_ref      = ref
)


"""
f,p = SA.get_power_spectrum(fs,SA.fft_plasma)
plt.figure(figsize=(10,5))
plt.plot(f,p)
plt.show()
"""

phi = SA.phi
phi2 = SA2.phi
phi3 = SA3.phi

"""
plt.figure(figsize=(10,5))
plt.plot(1000*time_tot*np.arange(0,len(phi))/len(phi),phi)
plt.show()
"""

#找到t_Ip在为-10和15的索引
begin = np.argmin(np.abs(t_Ip+10))
end = np.argmin(np.abs(t_Ip-20))

"""
plt.figure(figsize=(10,5))
plt.plot(t_Ip[begin:end],Ip[begin:end])
plt.show()
"""


#同时绘制Ip和phi,具有两个y轴,再画一条y=2pi的线
fig, ax1 = plt.subplots(figsize=(15,3))
ax1.plot(t_Ip[begin:end]+10,Ip[begin:end]/1e3,color='tab:red')

ax1.axvline(t_Ip[begin]+20,color='tab:gray',linestyle='--')
ax1.text(t_Ip[begin]+20.1,32,'Experiment Begin',color='tab:gray')
ax1.axvline(t_Ip[end]+5,color='tab:gray',linestyle='--')
ax1.text(t_Ip[end]+5.1,32,'Experiment End',color='tab:gray')


ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Ip (kA)', color='tab:red')
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.grid()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(1000*time_tot*np.arange(0,len(phi))[:int(0.25*len(phi))]/len(phi),density_string_integral(SA.phi_fix[:int(0.25*len(phi))]), color='tab:blue')
ax2.plot(1000*time_tot*np.arange(0,len(phi2))[:int(0.25*len(phi2))]/len(phi2),density_string_integral(SA2.phi_fix[:int(0.25*len(phi2))]), color='tab:orange')
ax2.plot(1000*time_tot*np.arange(0,len(phi3))[:int(0.25*len(phi3))]/len(phi3),density_string_integral(SA3.phi_fix[:int(0.25*len(phi3))]), color='tab:green')
ax2.legend(['sig1','sig2','sig3'])
ax2.set_ylabel('DensityStringIntegral/m^3', color='tab:blue')  # we already handled the x-label with ax1
ax2.tick_params(axis='y', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title('Shot#'+shot+':'+'Ip and DensityStringIntegral')
plt.savefig('Shot'+shot+'.png',dpi=300)
plt.show()

time_save = 1000*time_tot*np.arange(0,len(phi))[:int(0.25*len(phi))]/len(phi)
density = density_string_integral(SA.phi_fix[:int(0.25*len(phi))])
#将time_save和density保存到.txt文件中
np.savetxt('density.txt',np.c_[time_save,density])





"""
plt.figure()

plt.plot(SA.phi_fix)
plt.show()

#画出中频对应的ref和sig1
plt.figure(figsize=(10,5))
plt.plot(SA.fm_signal_plasma)
plt.plot(SA.fm_signal_ref)
plt.show()
"""