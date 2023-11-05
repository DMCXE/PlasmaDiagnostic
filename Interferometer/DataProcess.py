import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

class SignalAnalysis:
    def __init__(self,fs,fm_target,signal_plasma,signal_ref):
        self.fs             = fs                 #采样频率
        self.fm_target      = fm_target          #目标中频，考虑到实际情况可能存在频率移动
        self.signal_plasma  = signal_plasma      #等离子体信号
        self.signal_ref     = signal_ref         #参考信号
        self.lendata        = len(signal_plasma) #数据长度
        self.time_tot       = self.lendata/fs    #总时间

        #FFT
        self.fft_plasma = np.fft.fft(self.signal_plasma)
        self.fft_ref    = np.fft.fft(self.signal_ref)

        #提取信号
        self.fm_signal_plasma = self.get_fm_signal(self.fs,self.fm_target,self.fft_plasma)
        self.fm_signal_ref    = self.get_fm_signal(self.fs,self.fm_target,self.fft_ref)

        #鉴出频率
        self.phi = self.phase_detector(self.fs,self.fm_real,self.fm_signal_plasma,self.fm_signal_ref)
        self.phi_fix = self.phase_zero_fix(self.phi)

    def get_power_spectrum(self,fs,fft_signal):
        '''
        获得不同频率的归一化功率谱
        '''
        #通过fft获得结果以及频率分布
        frequencies = np.fft.fftfreq(len(fft_signal), 1/fs)

        #计算不同频率分量的功率谱密度
        power_spectrum = np.abs(fft_signal)**2

        #计算总功率
        total_power = np.sum(power_spectrum)

        #计算每个频率的占比
        power_percentages = (power_spectrum/total_power)*100

        #归一化
        power_percentages = power_percentages/np.max(power_percentages)

        return frequencies,power_percentages
    

    def get_main_frequencys(self,exist_frequencys=False,**kwargs):
        '''
        获得归一化功率谱中，占比最高的几个频率
        **kwargs:
            {fs:采样频率;signal:信号}
            或
            {frequencies:频率分布;power_percentages:功率谱密度}
        '''
        if exist_frequencys == False:
            fs = kwargs['fs']
            signal = kwargs['signal']
            frequencies,power_percentages = self.get_power_spectrum(fs,signal)
        else:
            frequencies = kwargs['frequencies']
            power_percentages = kwargs['power_percentages']
        
        main_frequencies = frequencies[power_percentages>0.05]
        return np.array(main_frequencies)
    
    def get_fm_real(self,fm_target,freq_main):

        #取出freq_main中，正的部分
        mains = freq_main[freq_main>0]
        fts = fm_target*np.ones(len(mains))
        index = np.argmin(np.abs(fts-mains))
        return mains[index]
    
    
    def get_fm_signal(self,fs,fm_target,fft_signal):
        
        freq,power_p = self.get_power_spectrum(fs,fft_signal)

        freq_main = self.get_main_frequencys(exist_frequencys  = True,
                                             frequencies       = freq,
                                             power_percentages = power_p)
        
        fm_real = self.get_fm_real(fm_target,freq_main)
        self.fm_real = fm_real

        freq_range = 0.3 * fm_real

        #确定对应频率范围内的索引
        indices = np.where(((freq >=   fm_real - freq_range) 
                         &  (freq <=   fm_real + freq_range))     #正频率部分
                         | ((freq <= -(fm_real - freq_range)) 
                         &  (freq >= -(fm_real + freq_range))))   #负频率部分   
        
        #滤波，只包含fm展宽内的频率分量，并由于DFT特点构造共轭
        filtered_fft = np.zeros(len(fft_signal),dtype=complex)
        filtered_fft[indices] = fft_signal[indices]

        fm_signal = np.fft.ifft(filtered_fft)

        return fm_signal

    def zero_detector(self,fs,signal):
        time_array = np.arange(0,self.time_tot,1/fs)
            #如果认为信号的相位差不会大于2pi，那么可以简化上述过程
        zero_crossing = np.where(np.diff(np.sign(signal))>0)[0]
        t = np.array([time_array[tt] for tt in zero_crossing])
        return t
        
    def phase_detector(self,fs,fm,signal_plasma,signal_ref):
        t_p = self.zero_detector(fs,signal_plasma)
        t_r = self.zero_detector(fs,signal_ref)
        return (t_r-t_p)*2*np.pi*fm
    
    def phase_zero_fix(self,phi):
        #对相位做零漂补偿
        phi_0 = phi[:int(0.1*len(phi))]
        phi_average = np.sum(phi_0)/len(phi_0)
        return phi - phi_average


def density_string_integral(phi,f=650*1e9):
    '''
    计算弦密度积分
    '''
    e = 1.60217662e-19
    m_e = 9.10938356e-31
    c = 299792458
    epsilon_0 = 8.854187817e-12
    return phi*np.pi*c*epsilon_0*m_e*f/(e**2)

       
if __name__ == '__main__':
    # 生成一个示例信号
    fs = 20000  # 采样频率为1 kHz
    fm_target = 30
    t = np.arange(0,30,1/fs)  # 1秒的时间序列
    freq1 = 10  # 10 Hz信号频率
    freq2 = 20  # 20 Hz信号频率
    signal = np.cos(2 * np.pi * freq1 * t) + 0.5 * np.cos(2 * np.pi * freq2 * t)
    def f(x):
        return np.cos(2 * np.pi * 10 * x + np.sin(x))  \
            + np.cos(2 * np.pi * 30 * x+np.sin(x))     \
            + np.cos(2 * np.pi * 1000 * x)             \
            + np.random.random(len(x))
    


    data = f(t)
    refdata = np.cos(2 * np.pi * 30 * t )


    sg = SignalAnalysis(fs,fm_target,data,refdata)
    
    fft_p = sg.fft_plasma
    fft_r = sg.fft_ref
    fm_p = sg.fm_signal_plasma
    fm_r = sg.fm_signal_ref
    t_p = sg.zero_detector(fs,fm_p)
    t_r = sg.zero_detector(fs,fm_r)

    plt.figure(figsize=(10,5))
    plt.plot(np.real(fft_p))
    plt.show()

    #绘制plasma和ref信号以及过零点
    plt.figure(figsize=(10,5))
    plt.plot(t,fm_p)
    plt.plot(t,fm_r)
    plt.plot(t_p,np.zeros(len(t_p)),'o')
    plt.plot(t_r,np.zeros(len(t_r)),'o')
    plt.show()


    plt.figure(figsize=(10,5))
    plt.plot(30*np.arange(0,len(sg.phi))/len(sg.phi),sg.phi,label='phi(t)',lw = 2)
    plt.plot(30*np.arange(0,len(sg.phi))/len(sg.phi),np.sin(30*np.arange(0,len(sg.phi))/len(sg.phi)),label='sin(t)',ls = '--',c = 'black',lw = 2)
    plt.grid()
    plt.legend()
    plt.show()


