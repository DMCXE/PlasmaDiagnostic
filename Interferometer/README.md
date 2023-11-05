# 等离子体诊断方法-第二章：折射

本文相关代码与介绍可在[PlasmaDiagnostic/Interferometer at main · DMCXE/PlasmaDiagnostic (github.com)](https://github.com/DMCXE/PlasmaDiagnostic/tree/main/Interferometer)中找到

## 目录
目录
1 多道偏振干涉仪诊断系统的应用
1.1 弦密度积分推导
1.2 探测器信号的采集与处理
1.2.1 探测器采样信号的构成
1.2.2 相位的提取
1.2.2.1 处理流程
1.2.2.2 离散傅里叶变换
1.2.2.3 中频提取
1.2.2.4 离散傅里叶逆变换
1.2.2.5 过零追踪鉴相器
1.2.2.6 零漂补偿
1.2.3 测试集测试
1.3 实验数据的处理结果
Shot#232091
Shot#232092
Shot#232094
1.4 误差分析与遗留问题
2. Homodyne & Heterodyne
3. 推导柱对称中的散射

# 1 多道偏振干涉仪诊断系统的应用

## 1.1 弦密度积分推导

要在冷等离子体假设下的折射系数的Appleton-Hartree公式：

$$
\mu^2=1-\frac{X(1-X)}{1-X-\frac12Y^2\sin^2\theta\pm[(\frac12Y^2\sin^2\theta)^2+(1-X)^2Y^2\cos^2\theta]^{1/2}}
$$

其中，  $X \equiv \omega_{p e}^2/\omega^2, \quad Y \equiv\omega_{c e}/\omega, \quad \mu \equiv k c/\omega$

在托卡马克诊断时，往往存在背景磁场，且波矢 $\vec{k}$ 与磁场  $B_0$ 相互垂直。对于寻常波入射，存在 $\vec{E}//\vec{B_0}$ ，此时 $Y=0$ , 这时，折射系数可以写为：

$$
\mu^2=1-X=1-\frac{\omega^2_{pe}}{\omega^2}=1-\frac{n_e}{n_c}
$$

光束经过等离子体产生的相移变换表示为：

$$
\phi_p = \int{(k_{plasma}-k_0})dl
$$

折射系数同样满足定义：

$$
\mu \equiv \frac{kc}{\omega}= \frac{k}{k_0}
$$

因此，相移可以进一步的化简为：

$$
\begin{aligned}
\phi_p = \int{(k_{plasma}-k_0})dl=-\int{\frac{\omega}{c}(1-\mu)}dl \\
=-\frac{\omega}{c}\int{(1-(1-\frac{n_e}{n_c})^{1/2})}dl
\end{aligned}
$$

考察$(1-n_e/n_c)^{1/2}$，当等离子体中电子密度 $n_e$  远远小于截止密度 $n_c$ 时，将此表达式应用泰勒展开有：

$$
(1-\frac{n_e}{n_c})^{1/2}=1-\frac12\frac{n_e}{n_c}-\frac18(\frac{n_e}{n_c})^2+o(...)
$$

截取到第二项，相移表达进一步简化为：

$$
\begin{aligned}
\phi_p =-\frac{\omega}{c}\int{(1-(1-\frac{n_e}{n_c})^{1/2})}dl=-\frac{\omega}{2n_cc}\int{n_e}dl
\end{aligned}
$$

截止密度$n_c$表示为(即对应截止频率时等离子体的振荡频率)：

$$
n_c=\frac{\varepsilon_0m_e}{e^2}\omega^2=\frac{\varepsilon_0m_e}{e^2}f^2[Hz]=\frac{\varepsilon_0m_ec^2}{e^2}\lambda^{-2}[m]
$$

移表达进一步简化为：

$$
\begin{aligned}
\phi_p =-\frac{\omega}{2n_cc}\int{n_e}dl=\frac{e^2}{\pi c\varepsilon_0m_e}f^{-1}\int{n_e}dl
\end{aligned}
$$

即得到了相位移与等离子体密度弦积分的关系。

## 1.2 探测器信号的采集与处理

### 1.2.1 探测器采样信号的构成

某干涉仪组成原理如下所示，是一个典型的外差干涉仪：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/46d98732-676d-49f1-9a08-701a021ef112)

在经过等离子体诊断区域的Beam1信号 $E_1$ 和LO信号 $E_2$ 可以表示为：

$$
\begin{aligned}
E_1 = A_1\cos(\omega_1t+\phi(t))\\
E_2 = A_2\cos(\omega_2t)
\end{aligned}
$$

未经过等离子体诊断区域的Beam1信号 $E_0$ 和LO信号 $E_2$ 可以表示为：

$$
\begin{aligned}
E_0 = A_0\cos(\omega_1t)\\
E_2 = A_2\cos(\omega_2t)
\end{aligned}
$$

在探测器上分别进行混频，有

$$
\begin{aligned}
\mathrm{signal\_plasma} =(E_1+E_2)^2 &=A_1^2\cos^2(\omega_1t+\phi)+A_2^2\cos^2(\omega_2t)\\
&+A_1A_2\cos[(\omega_1+\omega_2)t+\phi]+A_1A_2\cos[(\omega_1-\omega_2)t+\phi] \\
&=\frac12(A_1^2+A_2^2)+\frac12A_1^2\cos(2\omega_1t+2\phi)+\frac12A_2^2\cos(2\omega_2t)\\
&+A_1A_2\cos[(\omega_1+\omega_2)t+\phi]+A_1A_2\cos[(\omega_1-\omega_2)t+\phi]
\end{aligned}
$$

$$
\begin{aligned}
\mathrm{signal\_refer} =(E_0+E_2)^2 &=A_0^2\cos^2(\omega_1t)+A_2^2\cos^2(\omega_2t)\\
&+A_0A_2\cos[(\omega_1+\omega_2)t]+A_0A_2\cos[(\omega_1-\omega_2)t]\\
&=\frac12(A_1^2+A_2^2)+\frac12A_1^2\cos(2\omega_1t)+\frac12A_2^2\cos(2\omega_2t)\\
&+A_1A_2\cos[(\omega_1+\omega_2)t]+A_1A_2\cos[(\omega_1-\omega_2)t]
\end{aligned}
$$

对于此两束信号，可取出中频角速度为 $\omega_m=\omega_1-\omega_2$ 。在实际的诊断设备中，LO源与Beam源的频率可达到数百GHz，而两者的差频（中频）往往是数MHz，对探测器信号采样频率往往为数十MHz。

由**Nyquist采样定理**：

> 为了不失真地恢复模拟信号，采样频率应该不小于模拟信号频谱中最高频率的两倍。

因此在该混频信号中，可达到数百GHz的LO源与Beam源的频率在数十MHz的采样频率下，完全无法满足Nyquist定理，因而对应的信号全部失真。即原始混频信号中角速度为 $\omega_1,\omega_2,\omega_1+\omega_2$ 的的信号将失真，而信号差频（中频）对应的角速度分量 $\omega_m=\omega_1-\omega_2$ 由于只有数MHz，相对采样信号而言满足Nyquist定理，能够被不失真的被采集出来。

### 1.2.2 相位的提取

#### 1.2.2.1 处理流程

整体流程为

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/2b69a371-120d-4d2b-9148-c54795a5076b)

鉴相器组成为

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/976c46ff-d5c1-41b6-a90a-36f61417c770)


#### 1.2.2.2 离散傅里叶变换

对于实信号 $x(t)$ 的采样信号 $x[n]$ ，通过离散傅里叶变换DFT可以变换为 $X[t]$ ：

$$
X[k] = \sum_{n=0}^{N-1}x[n]e^{-j(2\pi/N)kn} \quad ,(0\le k \le N-1)
$$

对应的离散傅里叶逆变换为：

$$
x[n]=\frac1N \sum_{k=0}^{N-1}X[k]e^{j(2\pi/N)kn} \quad ,(0\le n \le N-1)
$$


需要注意的是：

1. 变换后序列 $X[k]$ 是复数序列，前一半数据与后一半数据共轭对称
2. 变换后序列 $X[k]$ 中k点的含义是： $f = k f_s/N$ ，即序列索引对应采样频率的 $1/N$

这也说明，如果需要通过离散傅里叶变换对频谱进行处理，最好在处理时同时考虑共轭频率部分对原始信号的贡献。

在python中， 这一步可以通过使用`np.fft`包中的`fft()`函数进行快速傅里叶变换（FFT）

```python
fft_signla = np.fft.fft(signal)
```

#### 1.2.2.3 中频提取

对原始信号进行带通滤波，在中频 $f_m$ 左右设置一个带通滤波器，提取中频 $f_m$ 范围内的主要频率，过滤掉高频与基频范围内的噪音。

在对于本数据的处理中，由于受到多普勒效应等或源信号频率误差，会导致实际中频并不严格的处于2.125Mhz，由于缺乏对实际物理图像的了解，导致中频会产生一点的偏差，在充分的了解频移特性与范围前，尚且无法很好的确定带通滤波器的范围。在此采用了一种伪带通滤波：先找到真实的中频，再依据一定展宽提取。

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/5e77897a-ab8c-4411-8e8a-900997248531)


在python中提供了多种便于计算处理DFT的工具：

根据`np.fft.fftfreq()`可以方便地获得每个索引对应的频率，由此即可计算归一化功率谱

```python
def get_power_spectrum(fs,fft_signal):
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
```

在归一化功率谱中，先选出占有主要功率强度的频率，通过与目标中频距离确定真实的中频位置：

```python
#确定主要功率强度
def get_main_frequencys():...
    return main_frequencies
#获得真实中频
def get_fm_real(fm_target,freq_main):...
	return fm_real
```

构建带通滤波器，获得中频分量

```python
freq_range = range * fm_real
return filtered_fft
```

#### 1.2.2.4 离散傅里叶逆变换

对于滤波后产生的频谱，进行离散傅里叶逆变换，即可获得中频信号

```python
fm_signal = np.fft.ifft(filtered_fft)
```

#### 1.2.2.5 过零追踪鉴相器

当对比信号和实际信号的中频被提取出来后，即可通过追踪在每个采样时间内正向通过零点的时间获得相位差。两信号正向过零的时间为：

$$
\begin{aligned}
\omega_m t_R = 2\pi m_R + \frac32\pi \\
\omega_m t_S+\phi_p = 2\pi m_S + \frac32\pi
\end{aligned}
$$

相移可以表示为：

$$
\phi_p = \omega_m(t_R-t_S)+2\pi(m_S-m_R)
$$

在粗略的认为相移不会超过$2\pi$的前提下，相移可以简写为

$$
\phi_p = \omega_m(t_R-t_S)
$$
正向过零点追踪能够被简写为：

```python
def zero_detector(fs,signal):
    time_array = np.arange(0,self.time_tot,1/fs)
    #如果认为信号的相位差不会大于2pi，那么可以简化上述过程
    zero_crossing = np.where(np.diff(np.sign(signal))>0)[0]
    t = np.array([time_array[tt] for tt in zero_crossing])
    return t
```

追踪出零点，即可给出在每个采样时间对应的相位差：

```python
def phase_detector(fs,fm,signal_plasma,signal_ref):
    t_p = self.zero_detector(fs,signal_plasma)
    t_r = self.zero_detector(fs,signal_ref)
    return (t_r-t_p)*2*np.pi*fm
```

#### 1.2.2.6 零漂补偿

由于诊断系统在空间中具有不对称性，信号采集器自身也可能存在温度部分不均等问题，导致诊断系统整体产生零点漂移，表现为在系统中为进行放电时，仍然能够测出相位差。因此需要对零点漂移进行修正。这里非常简单的考虑将相位差整体平移消除直流分量。

```python
 def phase_zero_fix(self,phi):
        #对相位做零漂补偿
        phi_0 = phi[:int(0.1*len(phi))]
        phi_average = np.sum(phi_0)/len(phi_0)
        return phi - phi_average
```

### 1.2.3 测试集测试

为了确定上述算法的可用性，这里通过一个示例：对于模拟数据

$$
f(t) = \cos(2\pi\cdot10\cdot t+\sin(t))+\cos(2\pi\cdot30\cdot t+sin(t))+\cos(2\pi\cdot1000\cdot t)+noise(t)
$$

对应参考频率信号为：

$$
f(t) = \cos(2\pi\cdot30\cdot t)
$$

采样频率为20Khz, 信号持续时间30秒，目标中频频率为30Hz。

经过上述流程，可以提取出中频信号与过零点：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/1488aaf6-9a38-4eba-a539-2ec4b615ab66)


并可以正确的提取出中频信号的相位差sin(t)：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/2bddabee-4bd1-4e0e-aa21-1a48ad3a3c9f)


验证了算法的正确性。

## 1.3 实验数据的处理结果

实验数据共有三组，分别为Shot#232091、Shot#232092和Shot#232094，每个数据文件中包含的信号种类如下图示

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/6f6cdfa4-e488-4f59-94bb-12508d832587)


数据处理时，首先指定文件(以sig1为例）：

```python
import scipy.io as sio
shot = '232094'
data = sio.loadmat(shot+'.mat')
Ip =data['ip1'].reshape(1,-1)[0]
t_Ip = data['t_ip1'].reshape(1,-1)[0]
sig1 = data['sig1'].reshape(1,-1)[0]
```

指定采集频率、目标中频

```python
fs        = 60   *1e6
fm_target = 2.125*1e6
```

将数据输入鉴相器

```python
SA = SignalAnalysis(
    fs              = fs,
    fm_target       = fm_target,
    signal_plasma   = sig1,
    signal_ref      = ref
)
```

将数据输入弦密度积分计算器

```python
phi = SA.phi_fix 					#获取零点漂移修正后的相位差数据
dsi = density_string_integral(phi)	#计算出弦密度积分
```

最后，能够得出Shot#232091、Shot#232091和Shot#232091的计算数据

### Shot#232091

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/b1f918e0-aaaa-4aba-8fbd-102bd90dc03f)

### Shot#232092

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/8220cab2-717d-41c2-9b91-a2c35921b11c)


### Shot#232094

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/df48a28c-f991-4389-9c50-bd61ddb4e105)


## 1.4 误差分析与遗留问题

1. 在消除零点漂移的时候，采用的方法过于简单，不能产生很好的效果
2. 实际上，注意到三炮数据在电流下降段，会出现多个密度峰值，在这些峰值处，实际对应的相位差大于 $2\pi$ ,实际上无法使用之前我们假设的相位差小于 $2\pi$ 的计算假设。但是如果归到 $[0,2\pi]$ 中，这些峰值又会变成瞬时的第密度甚至零密度。如何解释这些数据需要进一步讨论。
3. 本文数据处理的前提是：尽量模仿能够在FPGA或片上系统实现的信号处理方式。存在一些更先进的算法，例如Hilbert变换、EMD经验模态分解、小波变换等，具有很大的分析上述数据与解释的空间和准确性。为了更好的处理这些数据、获得更丰富的物理结果，希望能够在未来的工作中进一步发展。

# 2. Homodyne & Heterodyne

无论是零差干涉仪还是外差干涉仪，两臂信号混合之后均会产生两臂频率相加项和两臂频率相减相，主要关心频率相减项 $\omega_m = \omega_1 -\omega_2$ .

$$
\begin{aligned}
\mathrm{signal\_plasma} =(E_1+E_2)^2 &=A_1^2\cos^2(\omega_1t+\phi)+A_2^2\cos^2(\omega_2t)\\
&+A_1A_2\cos[(\omega_1+\omega_2)t+\phi]+A_1A_2\cos[(\omega_1-\omega_2)t+\phi] \\
&=\frac12(A_1^2+A_2^2)+\frac12A_1^2\cos(2\omega_1t+2\phi)+\frac12A_2^2\cos(2\omega_2t)\\
&+A_1A_2\cos[(\omega_1+\omega_2)t+\phi]+A_1A_2\cos[(\omega_1-\omega_2)t+\phi]
\end{aligned}
$$

$$
\begin{aligned}
\mathrm{signal\_refer} =(E_0+E_2)^2 &=A_0^2\cos^2(\omega_1t)+A_2^2\cos^2(\omega_2t)\\
&+A_0A_2\cos[(\omega_1+\omega_2)t]+A_0A_2\cos[(\omega_1-\omega_2)t]\\
&=\frac12(A_1^2+A_2^2)+\frac12A_1^2\cos(2\omega_1t)+\frac12A_2^2\cos(2\omega_2t)\\
&+A_1A_2\cos[(\omega_1+\omega_2)t]+A_1A_2\cos[(\omega_1-\omega_2)t]
\end{aligned}
$$

**零差探测（Homodyne  Detection）**

当 $\omega_m \ll \omega_1,\omega_2$ 时，即 $\omega_m = 0$ , 对应的探测技术即为零差探测技术。

经过等离子体的诊断信号对应有表达式：


$$
\begin{aligned}
\mathrm{signal\_plasma} =(E_1+E_2)^2 
&=\frac12(A_1^2+A_2^2)+\frac12A_1^2\cos(2\omega_1t+2\phi)+\frac12A_2^2\cos(2\omega_2t)\\
&+A_1A_2\cos[(\omega_1+\omega_2)t+\phi]+A_1A_2\cos(\phi)
\end{aligned}
$$

该信号的强度 $I$ 主要贡献项为:

$$
I_{plasma}\propto \frac12(A_1^2+2A_1A_2\cos(\phi)+A_2^2)
$$

对于对比信号

$$
\begin{aligned}
\mathrm{signal\_refer} =(E_0+E_2)^2 
&=\frac12(A_1^2+A_2^2)+\frac12A_1^2\cos(2\omega_1t)+\frac12A_2^2\cos(2\omega_2t)\\
&+A_1A_2\cos[(\omega_1+\omega_2)t]+A_1A_2\cos(0)
\end{aligned}
$$

该信号的强度$I$主要贡献项为:

$$
I_{refer}\propto \frac12(A_1^2+2A_1A_2+A_2^2)=\frac12(A_1+A_2)^2
$$

主要能够通过接受信号的强度判断出判断出相位变换。但是由于强度项具有平方特征，因此只能得到相位变化的绝对值，无法通过强度判断相位是正向变换还是反向变换，也就无法直接判断出等离子体弦密度积分是增长还是下降。

**外差探测（Heterodyne Detection）**

当 $\omega_m \not= 0$ , 对应的探测技术即为外差探测技术。当 $\omega_1,\omega_2$ 远远大于采样频率，而 $\omega_m$ 相比于 $\omega_1,\omega_2$ 较小，且相对于采样频率而言能够满足Nyquist采样定理（采样频率大于$f_m$两倍），此时就能够将混杂在信号中频率为 $f_m$ 的分量采集出来。即得到：

$$
\begin{aligned}
\mathrm{signal\_plasma} &=A_1A_2\cos[(\omega_1-\omega_2)t+\phi]\\
\mathrm{signal\_referen} &=A_1A_2\cos[(\omega_1-\omega_2)t]
\end{aligned}
$$

此时可以通过比较两个信号相位的相位差直接获得相位变化，且能够反应相位的相对变化，可以分辨出相位信号是正向的还是逆向的。即能够判断出等离子体弦密度积分是增长还是下降的。同时，外差探测不需要严格的保证两个源频率相同，只需要保证源本身输出频率是稳定的。同样的，如果只考虑相对相位的变化，通过合理的布置可以使其具有抗环境干扰能力，起到类似于差分补偿电路的作用。

# 3. 推导柱对称中的散射

**#存在大量问题，还需要进一步修改**

> Consider a beam propagating along a chord of a refractive cylinder at a distance y from the axis. Suppose the cylinder has a refractive index N(r), where r is the radius. Obtain a general equation for the angular deviation of the beam $\theta$ due to refraction when $\theta r\ll y$ so that the chord can be approximated as straight. In the case where $\omega\gg\omega_p$ and $n = n_o(1-r^2/a^2)$, calculate the value of y at which 0 is greatest and prove that this maximum 0 is $n_0/n_c$

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/d71990ce-1fc5-4581-9868-ca090bbd0994)


光线在介质中传播常见的方程形式为：

$$
\nabla n = \frac{d}{ds}(n\frac{dr}{ds})
$$

以圆柱圆截面圆心为原点建立极坐标系，由圆的的几何特性易得，射线方向与法线的夹角始终等于射线方向与径向方向的夹角。由射线折射的斯涅尔定律：

$$
n_1\sin(\theta_1) = n_2\sin(\theta_2)
$$

在球对称坐标系中，根据传播方程可以推出矢积关系：

$$
r \times \frac{dr}{d\eta} = const
$$

其中， $\eta$ 为路径L的弧长S ，该方程说明了光线始终在这个平面内传播，标量形式满足：

$$
N(r)\cdot r\cdot \sin(\alpha) = const
$$

其中，$\alpha$ 为射线方向与径向的夹角。在光线入射边界处，$N(R)=1$ , $R\sin(\alpha) = y$ ，则常数即为

$$
N(r)\cdot r\cdot \sin(\alpha) = y
$$

根据光线运动的几何关系，可以得到

$$
(\frac{r\cdot d\theta}{dr})^2 = \tan^2(\alpha)
$$

与折射定律联立，可以得到

$$
(\frac{r\cdot d\theta}{dr})^2 = \frac{1}{N^2(r)\cdot r^2/y^2-1}
$$

对偏转角进行积分, 其中 $r_0$ 是射线上距离原点最近的

$$
\theta(y) = \pi - 2\int_{r_0}^{\infin}\left|\frac{d\theta}{dr}\right|dr=\pi-2\int_{r_0}^{\infin}(N^2(r)r^4y^{-2}-r^2)^{-1/2}dr
$$

折射率与等离子体密度直接具有关系式：

$$
N(r) = (1-\frac{n_e}{n_c})^{1/2}
$$

当 $n_e = n_0(1-r^2/a^2)$ 时，有

$$
N(r) = \frac{r}{a}
$$

因此 $\theta(y)$ 的表达式进一步写为

$$
\begin{aligned}
\theta(y) = \sin^{-1}\frac{V_0\left[(\frac ya)^2-(\frac ya)^4\right]^{1/2}}{\left[V_0(\frac ya)^2+(\frac{1-V_0}{2})^2\right]^{1/2}}
\end{aligned}
$$

其中，  $V_0 = n_0/n_c $  。对上式求解，可以得到对于 $\theta$ 最大的值为 $V_0$


