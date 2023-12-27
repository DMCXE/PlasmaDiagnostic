import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import fsolve,brentq

def mcrossection(R0=1.85):
    # Inside wall
    pr = np.array([1.238, 1.238, 1.271, 1.399, 1.547, 1.718]) - R0
    pz = np.array([0, 0.916, 1.0447, 1.197, 1.254, 1.244])
    
    # Outside wall
    ppr = np.array([1.2, 1.2, 1.271, 1.399, 1.561, 1.742]) - R0
    ppz = np.array([0, 0.973, 1.14, 1.263, 1.316, 1.320])
    
    Qr = np.array([2.730, 2.84]) - R0

    # Limiter location
    nr = np.array([1.369, 1.369, 1.428, 1.328, 1.366, 1.594, 1.761, 1.822, 2.279, 2.65]) - R0
    nz = np.array([0, 0.46, 0.745, 1.016, 0.978, 1.103, 1.168, 0.926, 0.488, 0.488])
    
    plt.plot(ppr, ppz, linestyle='-', color='k', linewidth=2)
    plt.plot(ppr, -ppz, linestyle='-', color='k', linewidth=2)
    plt.plot(pr, pz, linestyle='-', color='k', linewidth=2)
    plt.plot(pr, -pz, linestyle='-', color='k', linewidth=2)
    
    plt.plot(nr, nz, linestyle='-', color='b', linewidth=2)
    plt.plot(nr, -nz, linestyle='-', color='b', linewidth=2)

    # Additional plots
    ap = Qr[0] - pr[-1]
    bp = pz[-1]
    app = Qr[1] - ppr[-1]
    bpp = ppz[-1]
    theta = np.linspace(np.pi/2, 0, 50)
    
    pR = pr[-1] + ap * np.sin(theta)
    pZ = bp * np.cos(theta)
    
    plt.plot(pR, pZ, linestyle='-', color='k', linewidth=2)
    plt.plot(pR, -pZ, linestyle='-', color='k', linewidth=2)
    
    ppR = ppr[-1] + app * np.sin(theta)
    ppZ = bpp * np.cos(theta)
    
    plt.plot(ppR, ppZ, linestyle='-', color='k', linewidth=2)
    plt.plot(ppR, -ppZ, linestyle='-', color='k', linewidth=2)

    # 修正这一行的绘制方式
    plt.plot(np.zeros(100), np.linspace(-max(ppz), max(ppz), 100), linestyle='--', color='black', linewidth=1)

    plt.xlabel('Radial Location  /m', fontsize=14, fontname='times new roman')
    plt.ylabel('Vertical Location  /m', fontsize=14, fontname='times new roman')

    plt.axis([min(ppr)-0.05, max(pR)+0.05, min(-ppz)-0.1, max(ppz)+0.1])

    plt.axis('equal')
    p1 = min(-ppz)-0.1
    p0 = 0+0.2
    plt.text(p0, p1, 'Magnetic Axis', fontsize=12, verticalalignment='bottom', horizontalalignment='right')

def plasma_cross(R0,a,b,kappa,delta,Npol = 100):
    theta = np.linspace(0,2*np.pi,Npol)
    R = R0 - b + (a + b*np.cos(theta))*np.cos(theta + delta*np.sin(theta))
    Z = kappa*a*np.sin(theta)
    plt.plot(R,Z,'--',c='purple')




if __name__ == '__main__':
    e  = 1.602e-19
    me = 9.109e-31

    plt.figure(figsize=(6,12))
    mcrossection()
    
    R0    = 1.85
    B0    = 2.2486
    n     = 2
    gamma = 1
    B_func      = lambda R   : B0*R0/R
    omega_func  = lambda B   : n*e*B/me

    omega_setting = np.arange(105.1,167.1+2,2) * 1e9 #f: GHZ
    observation_points = np.zeros_like(omega_setting)

    # 通过fsolve进行正向解，ω(R)-2πf_target = 0
    for index,omega in enumerate(omega_setting):
        #使用fsolve方法
        R_func = lambda R : omega_func(B_func(R)) - omega*2*np.pi
        observation_points[index] = fsolve(R_func,1.4)

    print(observation_points)
    

    # 通过公式进行反向解
    observation_points_func = lambda f: n*B0*R0*e/(2*np.pi*me*f)
    observation_points = observation_points_func(omega_setting)
    print(observation_points)

    # 由于R0在截面图中的位置为0，计算的观察点需要对R0修正
    opo = observation_points - R0
    plt.scatter(opo,np.zeros_like(omega_setting),c='red',s=8)

    plasma_cross(0,R0-0.999*observation_points[-1], 0.28,1.4,0.0)
    plasma_cross(0,R0-0.99*observation_points[-5], 0.20,1.4,0)
    plasma_cross(0,R0-0.99*observation_points[-10], 0.15,1.4,0)
    plasma_cross(0,R0-0.99*observation_points[-15], 0.05,1.4,0)
    plasma_cross(0,R0-0.99*observation_points[-25], 0,1.4,0)

    plt.savefig("cross.png",dpi=300,bbox_inches = 'tight') 
    plt.show()


[2.21576223, 2.17438479, 2.13452439, 2.09609911, 2.05903281, 2.02325465,
 1.98869864, 1.9553032 , 1.92301082, 1.89176776, 1.86152367, 1.8322314,
 1.80384671, 1.77632808, 1.74963644, 1.72373509, 1.69858943, 1.67416686,
 1.65043665, 1.62736975, 1.60493874, 1.58311768, 1.56188203, 1.54120854,
 1.52107518, 1.50146106, 1.48234634, 1.4637122 , 1.44554072, 1.4278149,
 1.41051854, 1.39363621]

# 题目三