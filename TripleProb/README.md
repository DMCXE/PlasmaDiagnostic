# 等离子体诊断方法-第一章：静电探针
# 写在前面——代码基本用法
## 库需求

scipy—用于.mat文件读取

numpy—基础运算库

matplotlib—绘图与动画生成

tqdm—绘制进度条

可选：

ffmepg—生成mp4的工具库

multiprocessing—python的并行库，用于数据的并行处理与保存

如果需要可通过

```bash
pip install scipy numpy matplotlib tqdm multiprocessing
```

或者通过conda环境

```bash
conda install scipy numpy matplotlib tqdm multiprocessing
```



## 基本用法

### 数据处理并生成处理后文件

在文件 `data_process.py` 中的 `if __name__ == '__main__':`入口除，对于函数`per_processing_all_data(mutiprocess = True)`，需要手动将函数中文件路径修改为本机上原始数据的代码：

```python
def per_processing_all_data(mutiprocess = False):
	...
        for d1 in np.arange(1,14):
            for d2 in np.arange(1,26):
                data_path = './dataall/dataAll/y'+str(d1)+'x'+str(d2)+'/'
                path_list.append(data_path)
       ...
```

将 `data_path = './dataall/dataAll/y'`改为储存中的全部数据的路径。

之后，在函数`def position_average`中，将最后`scio.savemat('./dataall/data_processed/'+...)`路径修改为需要保存的位置。最后，在入口处：

```python
if __name__ == '__main__':
    st = time.time()
    per_processing_all_data(mutiprocess = True)
    stt = time.time()
    print(stt-st)
```

如果不需要并行运算，mutiprocess项可选为False。

### 绘制处理后文件

在文件`data_ansys.py`中，首先在`class DrawAndAnimation()`中，找到函数` def Tensor_generate`, 并将其中data_path修改为上一节中处理后文件的保存路径。

如果需要绘制包含电子温度、电子密度与磁场强度与磁力线分布的动画并保存，在程序入口处执行：

```python
if __name__ == "__main__":
   DA = DrawAndAnimation()
   DA.BigAnimation(ways = 'x1y4')
```

其中，ways表示绘制出图像的排布，表示xy方向的图片排布数量。

如果只需要绘制电子温度与电子密度的童话，执行

```python
if __name__ == "__main__":
   DrawTwoatOnce()
```

这一步必须要求multiprocessing。


## 目录

等离子体诊断方法-第一章：静电探针 <br />
目录 <br />
1 三探针的推导 <br />
--1.1 三探针数据的推导 <br />
--1.2 三探针中解的存在性 <br />
2 磁探针的问题 <br />
--2.1 磁探针的标定 <br />
--2.2 双绞线的好处 <br />
3 LMP中的数据处理 <br />
--3.1 通过探针数据获得物理量的原理 <br />
----3.1.1 电子温度与电子密度数据的获得与计算 <br />
----3.1.2 磁场数据的处理 <br />
--3.2 数据的预处理 <br />
----3.2.1 驱动电流的相位误差 <br />
----3.2.2 创建驱动电流时间窗 <br />
------3.2.2.1 寻找驱动电流起始点与结束点的方法 <br />
------3.2.2.2 创建诊断时间窗 <br />
----3.2.3 单独位置数据的处理 <br />
----3.2.4 数据的存在性 <br />
------3.2.4.1 三探针诊断中的数据 <br />
------3.2.4.2 磁探针诊断中的数据 <br />
--3.3 数据预处理的程序介绍 <br />
----3.3.1 程序流图 <br />
----3.3.2 数据结构 <br />
--3.4 磁重联二维参数的演化 <br />
4 思考题 <br />


# 1 三探针的推导

## 1.1 三探针数据的推导

为了广泛的讨论，并依据文献[Chen, S., & Sekiguchi, T. (1965). Instantaneous Direct‐Display System of Plasma Parameters by Means of Triple Probe. Journal of Applied Physics, 36(8), 2363–2375. https://doi.org/10.1063/1.1714492]的相关假设，将三探针一般化为带有对称偏压的探针电路，后续将通过分析该电路以判断三探针诊断的适用条件。

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/45ea4483-f23c-4925-aa66-13b805b4993e)

由基尔霍夫定律，该电路系统的电压与电流方程可以写为：

$$
I_1=I_2+I_3
$$

$$
V_{d2}=V_2-V_1
$$

$$
V_{d3}=V_3-V_1
$$

对于单独的静电探针，以流入探针的电流（离子流）为正方向，由于每个静电探针的特性是相同的，均由电子流与离子流构成。其中，离子流是收到探针电位驱动的离子饱和流，电子流是麦克斯韦分布的热扩散流。电子流可以写为：

$$
J_e = en_e\sqrt{\frac{kT_e}{2\pi m_e}}
$$

因此此三个探针的电流方程可以写为：

$$
-I_1 = -SJ_{e}\mathrm{exp}(-\frac{eV_1}{kT_e})+SJ_i(V_1)
$$

$$
I_2 = -SJ_{e}\mathrm{exp}(-\frac{eV_2}{kT_e})+SJ_i(V_2)
$$

$$
I_3 = -SJ_{e}\mathrm{exp}(-\frac{eV_3}{kT_e})+SJ_i(V_3)
$$

其中，S是探针的表面面积，这里也假设是相同的。

假设由于电位变化引起的离子饱和流 $J_i(V)$ 变化远远小于电子电流的变化，即可认为 $J_i(V_1)\approx J_i(V_2)\approx J_i(V_3)$ ，对于上述电流方程，相加消去 $J_i(V)$ 后可以得到：

$$
\frac{I_1+I_2}{I_1+I_3}=\frac{1-\mathrm{exp}(-e(V_2-V_1)/kT_e)}{1-\mathrm{exp}(-e(V_3-V_1)/kT_e)}=\frac{1-\mathrm{exp}(-eV_{d2}/kT_e)}{1-\mathrm{exp}(-eV_{d3}/kT_e)}
$$

当采用PPT中所使用的简化三探针后：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/be1e40f5-8efe-496a-a004-50b5b0e06144)


此时，可认为$V_{d_2}$断路，便有 $I_2 = 0$,$I_1=I_3$ ,  $V_{d2}=V_f-V_1$ ，$V_{d}=V_3-V_1$ ，电流比值的表达式即可写为：

$$
\frac{1-\mathrm{exp}(-e(V_f-V_1)/kT_e)}{1-\mathrm{exp}(-eV_d/kT_e)}=\frac{I_1+0}{I_1+I_3}=\frac12
$$

当 $V_d \gg kT_e/e$ 时， $\mathrm{exp}(-V_d/kT_e)\approx 0$ 上式可以简化为:

$$
\mathrm{exp}(-e(V_f-V_1)/kT_e) = \frac 12
$$

化简后即为：

$$
T_e=\frac{e}{k\mathrm{ln}2}(V_1-V_f)
$$

对于电子密度，在探针的鞘层边界上将等离子体区域划分为鞘层区域、准中性等离子体区域与等离子体区域，如图示：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/29b590e1-c780-48ad-b738-90d5826eb14b)

理论上，准中性等离子体区域对应的电势为

$$
V_0 = \frac{kT_e}{2e}
$$

在鞘层与准中性等离子体区域边界上,有 $T_e\gg T_i$ ,离子速度与饱和离子流密度可以表示为:

$$
v_{is}=(2eV_0/m_i)^{1/2}=(kT_e/m_i)^{1/2}
$$

$$
J_{is}=en_{is}v_{is}=en_{is}(kT_e/m_i)^{1/2}
$$

在鞘层外的准中性区中，可以认为离子密度与电子密度近似相同，且电子密度分布符合麦克斯韦分布

$$
n_{is}\approx n_{ne}= n_e\mathrm{exp}[-\frac{e}{kT_e}\cdot \frac{kT_e}{2e}]=n_e\mathrm{exp}(-\frac12)\approx 0.61n_e
$$

则离子饱和流密度表示为：

$$
J_{is}=en_{is}v_{is}=0.61en_e(kT_e/m_i)^{1/2}
$$

对于面积为 $A_p$ 的探针，离子饱和流为：
$$
I_i = 0.61en_eA_p(kT_e/m_i)^{1/2}
$$

## 1.2 三探针中解的存在性

对于一般三探针的表达式：

$$
\frac{I_1+I_2}{I_1+I_3}=\frac{1-\mathrm{exp}(-e(V_2-V_1)/kT_e)}{1-\mathrm{exp}(-e(V_3-V_1)/kT_e)}=\frac{1-\mathrm{exp}(-eV_{d2}/kT_e)}{1-\mathrm{exp}(-eV_{d3}/kT_e)}=const
$$

令 $\phi=e/kT_e$ ，求解电子温度 $T_e$ 等价于求解 $\phi$ ，而 $\phi$ 的存在性可以表示为函数 $f(\phi)$ ：

$$
f(\phi)=c\mathrm{e}^{-\phi V_{d3}}-\mathrm{e}^{-\phi V_{d2}}+1-c=0

$$
存在非零解。注意观察易知函数 $f(\phi)$ 必然存在解 $\phi=0$ 。因此函数 $f(\phi)$ 存在非零解等价于存在极值点（这是一个必要不充分条件）

$$
f'(\phi)=-cV_{d3}\mathrm{e}^{-\phi V_{d3}}+V_{d2}\mathrm{e}^{-\phi V_{d2}}=0
$$

若极值点，则解得

$$
\phi = \frac{\mathrm{ln}(c\cdot V_{d3}/V_{d2})}{V_{d3}-V_{d2}}
$$

解的存在性要求了

$$
V_{d3}\cdot V_{d2}>0
$$

该表达式要求了非公共电极上必须具有与公共电极相比相同方向的电势差。在简化的三探针模型中，由于 $V_2$ 电位为低电位，即要求了悬浮电位$V_f$同为低电位，必须低于正电位 $V_1$ 。在后文的数据处理中，我们会遇到违反这一要求的情况。在3.4.2.1中我们将着重讨论这个问题。

# 2 磁探针的问题

## 2.1 磁探针的标定

磁探针的感应电动势为：

$$
\mathscr{E} = -\frac{d\phi}{dt} = -S_{eff}\frac{dB_{\bot}}{dt}
$$

其中， $S_{eff}$ 表示为线圈的有效截面积

对于探针有效截面 $S_{eff}$ 难以直接测量，需要进行标定，手段有：

1. 二次标定：将待标定探针与尺寸已知的标准线圈在同一个脉冲磁场中比较输出。标准线圈的尺寸不必很小。
2. 直接标定：将待标定的探针放置于强度已知的脉冲磁场中，测量输出电压。这种方法可以直接标定整个探针系统。
3. 直接测量：直接测量线圈的几何尺寸，包括直接测量线圈本身与引线等面积。

标定用的标准磁场可以用通过单层长螺线管线圈与亥姆霍兹线圈的脉冲电流I产生。

## 2.2 双绞线的好处

磁探针的结构如下图所示：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/06f9abb9-4a3b-4020-8812-f46fb10b00b6)

双绞线是一对相互绝缘的导线，按照一定规律相互缠绕绞合。由于磁探针的数据来源于等离子体系统中变化的磁通量，在传输线中同样会产生不同频率的电流。如果使用普通单线传输，这些电流会在系统中产生电磁干扰，产生电缆间的串扰，增加数据的噪声。使用双绞线是，由于绞线中同时传递的是幅值相等极性相反的信号，传输线路中的噪声会产生一组共模信号，一条导线上产生的干扰会被另一根上的干扰抵消。在等离子体诊断中的高能量、强磁场、多扰动的情况下可以起到很好的抗干扰作用，提升信号结构的完整性，提升不失真前提下采集信号的频率。同时，双绞线的使用与铺设成本还比较低，在诊断这种极端环境中，易于更换。

更优化的选择，一方面可以考虑在已有双绞线的基础上铺设金属网屏蔽层，进一步的加强抗电磁干扰能力。另一方面可以更换为同轴电缆，其具有多种绝缘包层，能够进一步的提升抗干扰能力。

部分观点参考自：高温等离子体诊断技术, 上海科学技术出版社

# 3 LMP中的数据处理

在本文中，获取数据的策略是使用单组三探针设备，在步进电机控制下，对y方向上13组，x方向上25组，共13*25个点进行诊断，每个点均诊断10次以尝试消除误差，每次诊断时其余参数均保持不变。数据结构为：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/f873b005-3014-42ff-b90f-61bcb4eac02f)

## 3.1 通过探针数据获得物理量的原理

### 3.1.1 电子温度与电子密度数据的获得与计算

根据三探针的相关计算，获得的电子温度通过上述数据可以表示为：

$$
T_e = \frac{e}{k\mathrm{ln}2}(V_+-V_f)
$$

电子密度可以通过离子饱和流 $I_{is}$ 获得：

$$
I_{is}=0.61en_eA_p\sqrt{kT_e/m_i}
$$

为了简化数据处理过程，结合数据结构中的离子饱和流信号（复合） $V_{ion}$ ，这里将电子密度表示为：

$$
n_e = 0.909\cdot 10^{19}\cdot V_{ion}
$$

### 3.1.2 磁场数据的处理

磁场诊断采用的是内部磁探针。内部磁探针实际上是一个插入等离子体的小螺线管线圈，利用了变化磁通量产生的感生电动势：

$$
\mathscr{E}=-\frac{d\phi}{dt}=-S_{eff}\frac{dB_{\bot}}{dt}
$$

其中， $S_{eff}$ 为标定的有效面积。根据实验中的相关规范，标定的二维磁探针有效面积为：x方向: 3.001691131e-4, y方向：4.550113814e-4。

在本实验的数据处理中，通过对xy方向的磁探针信号进行积分既可以得到。需要注意的是，通过磁探针的计算表达式不难发现，仅当磁场处于非稳态时，磁探针数据才可以通过积分获得磁场的相对变化量。除非我们知道磁场的初始值，否则只能分析它的涨落情况。

## 3.2 数据的预处理

### 3.2.1 驱动电流的相位误差

由于设备本身存在的固有精度等问题，即使设定上保证了每次放电参数、时间均一致，但实际诊断数据中，在一定精度上任然能发现误差。在2.a.iii中，规定了重联发生的时序信号以重联驱动电流为准，并取三段式放电中电流增长起始点作为重连开始时刻。但由于上述精度问题，实际诊断数据中可以发现驱动电流起始点时序并不完全一致，以data文件中，y4x15, y4x16, y7x14, y7x15, y9x9, y9x10, 此六个位置共60个数据的磁重联驱动电流为例，如图示：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/de4b1de5-98fe-4d48-af66-34a57401b128)

驱动电流幅值上误差较小，相位上误差较大。驱动电流的起始到结束过程分别对应磁重联的启动与结束，因此为了合理的评估并绘制磁重联过程中的诊断因素，需要依据驱动电流对齐时间序列，再对数据进行进行平均等操作。

### 3.2.2 创建驱动电流时间窗

LMP单次放电诊断全过程可以大致分为三个阶段：等离子体建立-平衡阶段、磁重联阶段、磁重联结束重新建立平衡阶段。以重联驱动电流发生的时序信号为基准时间，可以通过对每组数据创建时间滑窗的方式消除相位差。由于本实验针对的主要物理过程是磁重联过程，因此对于前后阶段的考量可以忽略，并且时间由绝对时间转换至磁重联阶段的相对时间。主要的流场可以概括为：

1. 计算一个采样点内10次的诊断数据
2. 创建时间窗口，仅包含重连发生窗口区间
3. 对齐重连起始时间
4. 平均参数值，将重连起始处即为0ms

#### 3.2.2.1 寻找驱动电流起始点与结束点的方法

即使是再驱动电流产生前的零电流阶段，由于背景噪声等原因，数据并不是平滑单调的曲线，如图所示：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/e61fee45-726a-4212-964a-aad77ef4a421)

这为采用一些成熟的寻找拐点的数学工具制造了一些困难，往往无法准确的找到电流上升的起始点与结束点。但是，由于采用的等时间步长的，因此在电流起始点后，单时间步长对于的电流增长远大于噪音涨落，即该数据的一阶差分的幅值能够帮助我们找到拐点。以y1x12/N1的驱动电流差分值为例：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/b1939f94-d36e-4291-9362-4612a32a1a87)

通过指定一阶差分中阈值，即可找到驱动电流起始点与结束点，并且时间相位上的误差仅存在于几个时间步长上。通过简单的构造，能够得到上述电流簇的每一个电流起始点的位置。当需要获得驱动电流结束点时，仅需要将电流序列反向后即可。

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/54d503f6-e3c6-4c74-8ade-f487fb9baccd)

```python
def find_kneed_arg(self, data1d):
        """
        找到上升拐点的位置
        data1d: 一维数据
        s: 为阈值
        return 拐点的位置索引
        """
        diff = np.diff(data1d)
        point = np.where(diff > self.s)
        point = point[0][0]
        return point
```

#### 3.2.2.2 创建诊断时间窗

无论各位置点驱动电流起始点在那个时刻，当规定起始点与结束点后，即可创建任意长度的时间窗，并按照此时间窗截取对于组数据中其它我们感兴趣的参数。这样即可构建出任意组具有完全相同长度、完全相同磁重联起始点与结束点的诊断数据。

创建时间窗后，数据如下所示：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/7bcc2bd7-d752-45e8-9d2c-b830423d21f4)

```python
def snap_data(self, data1d, kneed_point, before_kneed, after_kneed):
        """
        将各个数据依据磁重联起始点与过程创建时间窗以对齐时间测量误差
        data1d          :   一维数据
        kneed_point     :   磁重联脉冲起始点
        before_kneed    :   磁重联起始点前选取的步长
        after_kneed     :   磁重联起始点后选取的步长
        """
        data_window1d = data1d[kneed_point - before_kneed:kneed_point + after_kneed]
        return data_window1d
```

### 3.2.3 单独位置数据的处理

当创建时间窗后，才可对数据进行平均。若是直接对数据平均，收到相位差的影响会使得参数整体程序较低的趋势，采用直接平均与窗口后平均产生的数据差别如下图所示：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/ce56a671-ecaf-41dc-b703-2edd3ed6166b)

因此对于参数，需要先创建时间窗后再平均数据，以确保数据本身的固有参数特征被保留。在此流程下，即可对单独位置重复采集数据的平均。

```python
    def position_average_datas(self,datas):
        """
        对N维位置数据进行平均，采用np.mean的形式，可拓展
        datas: 高维位置数据
        """
        data_average = np.mean(datas)
        return data_average
```

### 3.2.4 数据的存在性

#### 3.2.4.1 三探针诊断中的数据

在1.2中，我们讨论过三探针中解的存在性。根据电子温度与电子密度的表达式：

$$
T_e = \frac{e}{k\mathrm{ln}2}(V_+-V_f)
$$

$$
I_{is}=0.61en_eA_p\sqrt{kT_e/m_i}
$$

可以发现，当 $V_+ < V_f$ 时， $T_e < 0$ ，此时温度与密度均无意义。实际实验中，往往通过对三探针施加很大的正电位以使得 $V_+ > V_f$ 始终成立，但分析本实验的数据，可以发现在磁重联启动期间，存在 $V_+ < V_f$ ，以y1x12/N1数据为例：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/dee618da-cd1b-4826-820b-dd2ef02308e3)

在磁重联建立后，出现了 $V_+ < V_f$ 的情况，此时对应的温度与密度为：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/b2428a8f-08be-45d4-b193-94d11f9e141c)

显然，当 $V_+ < V_f$ 关系被违反后，电子温度项与密度项均出现了负值。但是我们没办法处理这件事，甚至最多只能对温度项进行处理。在LMP的规范中，电子密度被写为了一个笼统的表达式：

$$
n_e = 0.909\cdot 10^{19}\cdot V_{ion}
$$

当 $V_+ < V_f$ 时，我们可以简单的认为此时**受到到磁重联的磁场配置，使得诊断区域的等离子体被排开，密度的大量降低导致探针局部不满足准中性条件，探针中收集的电流与电位信息不满足推出上述公式过程中基本假设**

在后续的数据分析中，我们允许出现负值情况，以对全局进行定性考量。 

#### 3.2.4.2  磁探针诊断中的数据

磁探针诊断数据中，需要明确的物理过程是：主要的磁场由磁重联驱动电流产生，即在磁重联驱动电流产生前后，磁探针测得的电势数据来源于等离子体的建立与涨落，不具有较大的物理意义。因此在后续分析中，采取手动归零的方式将磁重联驱动电流时间窗外的全部磁探针数据归零。 

## 3.3 数据预处理的程序介绍

### 3.3.1 程序流图



![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/8cedb823-45af-4b5e-9bf7-133877cde2ef)

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/b15cbf33-5ac8-49f6-9c6e-9943bee39965)

### 3.3.2 数据结构

数据的预处理由以下代码构成：

直接获取诊断数据信息并进行数据计算

```python
class Magnetic_Reconnection_Diagnostic :
    def __init__(self,data_path,B_Diagnostic = True,)
    def electron_temperature(self):...
 	def electron_density(self):...
	def magnetic_field_integral(self,B_data_1d,delta_time,s_eff,B0 = 0):...
	def magnetic_signal_line_average(self,B_data_1,B_data_2):...
```

获得数据拐点并生成时间窗：

```python
class DataProcess:
    def __init__(self,data_origin,find_kneed_error = 0.006,window_step = [100,500]):
    def find_kneed_arg(self, data1d):...
    	"""找到上升拐点的位置"""
    def find_kneed_down_args(self,data1d):...
    	"""找到下降拐点的位置"""
    def snap_data(self, data1d, kneed_point, before_kneed, after_kneed):...
    	"""将各个数据依据磁重联起始点与过程创建时间窗以对齐时间测量误差"""
    def position_average_datas(self,datas):
        """对N维位置数据进行平均，采用np.mean的形式，可拓展"""
```

对单独位置进行多次测量数据的平均

```python
def position_average(mkr_path = './data/y4x15/'):
    MRD = Magnetic_Reconnection_Diagnostic(data_path,
                                               B_Diagnostic=True)
    DP = DataProcess(MRD.I_drive,window_step=[100,700])
    Te_window = MRD.Te[window_index[0]:window_index[1]]
    ne_window = MRD.ne[window_index[0]:window_index[1]]
    B_x_m_ig = MRD.magnetic_field_integral(Bx,MRD.delta_t,s_eff=3.001691131e-4)
    B_y_m_ig = MRD.magnetic_field_integral(By,MRD.delta_t,s_eff=4.550113814e-4)
    average()...
    savedata...
```

通过并行化生成每一个诊断点上的数据。

```python
def per_processing_all_data(mutiprocess = False):
    if mutiprocess == True:
        import multiprocessing as mp
        from multiprocessing import RLock
        '''生成包含d1,d2的path的列表'''
        path_list = []
        for d1 in np.arange(1,14):
            for d2 in np.arange(1,26):
                data_path = './dataall/dataAll/y'+str(d1)+'x'+str(d2)+'/'
                path_list.append(data_path)
        '''多进程处理'''
        with mp.Pool(processes=12,
                     initializer=tqdm.tqdm.set_lock,
                     initargs=(RLock(),)
                     ) as pool:
            M = list(tqdm.tqdm(pool.imap(position_average,path_list),total=len(path_list)))
    else:
        for d1 in np.arange(1,14):
            for d2 in np.arange(1,26):
                data_path = './dataall/dataAll/y'+str(d1)+'x'+str(d2)+'/'
                position_average(data_path)
```

处理后每个位置生成的新文件储存的字典与规范为：

```python
'Te_average':Te_average,
'ne_average':ne_average,
'I_driven_average':I_driven_average,
'Bxmi_average':Bxmi_average,
'Bymi_average':Bymi_average
```

## 3.4 磁重联二维参数的演化

磁重联过程的动画将在附录文件里面展示，这里展示重要部分的静态图片。在附录的动态文件中，将详细的给出电子密度、电子温度、磁场分布以及磁力线分布随电流演化的情况。

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/1a3d241b-87a4-4355-8425-2bf44a947e66)


在静态演化展示中， 我们给出了电子温度与电子密度演化的大致趋势，并浮上来与磁力线相关的分布变化的情况。

对于电子密度随时间演化：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/4a284f92-951c-4d97-a6f1-0f6f90b51cd2)

对于电子温度随时间演化：

![image](https://github.com/DMCXE/PlasmaDiagnostic/assets/61933011/b86b345a-5d6b-4910-8781-1adf2d97e19e)


# 4 思考题

在数据分析中，我们已经得到了实验数据的构成：每一个位置测量十遍，在截面上测量25*13个位置，即6250个数据，每个数据唯一对应一次实验，即不可能在同一时刻对所有位置进行诊断。这样处理数据的基础在于等离子体的某些性质在相对较短的时间内（相对于等离子体的时间尺度）保持稳定，比如说从一个状态到另一个稳定状态时的弛豫时间，探针的响应要足够的快，以确保不会忽略重要的细节。这允许将一系列在不同时间点测量的数据视为同一个物理过程的时间演化。

在每次诊断时，需要确保在同一位置不同时间点进行测量时，等离子体的温度、密度、磁场等保持相对稳定，以减小测量时变化。对于探针，需要保证在每次诊断时，对于探针的标定与校准是一致的。
