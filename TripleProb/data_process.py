import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import tqdm
import time


class Magnetic_Reconnection_Diagnostic :
    def __init__(self,
                 data_path,
                 B_Diagnostic = True,):
        f                   = 2e6 # 2MHz
        self.data_path      = data_path
        self.data           = scio.loadmat(data_path)['data']
        self.time_step      = len(self.data)
        self.time           = 1e3 * np.arange(0, self.time_step, 1) / f   #ms
        self.delta_t        = 1 / f #ms 采样时间间隔

        if B_Diagnostic == True:
            self.B_x1       = self.data[:, 3]       #x方向双绞线1磁场信号
            self.B_x2       = self.data[:, 4]       #x方向双绞线2磁场信号
            self.B_y1       = self.data[:, 5]       #y方向双绞线1磁场信号
            self.B_y2       = self.data[:, 6]       #y方向双绞线2磁场信号

            #self.B_x1_ig = self.magnetic_field_integral(self.B_x1,self.delta_t,s_eff=3.001691131e-4)
            #self.B_x2_ig = self.magnetic_field_integral(self.B_x2,self.delta_t,s_eff=3.001691131e-4)
            #self.B_y1_ig = self.magnetic_field_integral(self.B_y1,self.delta_t,s_eff=4.550113814e-4)
            #self.B_y2_ig = self.magnetic_field_integral(self.B_y1,self.delta_t,s_eff=4.550113814e-4)

            self.B_x_mean = self.magnetic_signal_line_average(self.B_x1,self.B_x2)
            self.B_y_mean = self.magnetic_signal_line_average(self.B_y1,self.B_y2)

            self.B_x_m_ig = self.magnetic_field_integral(self.B_x_mean,self.delta_t,s_eff=3.001691131e-4)
            self.B_y_m_ig = self.magnetic_field_integral(self.B_y_mean,self.delta_t,s_eff=4.550113814e-4)


        self.V_source       = self.data[:, 0]       #阴极源放电电压  V
        self.I_cathode      = self.data[:, 1]       #阴极放电电流   I,kA
        self.I_drive        = self.data[:, 2]       #重联驱动电流   I,kA
        self.V_ion          = self.data[:, 7]       #离子饱和流信号         V_ION
        self.V_plus         = self.data[:, 8] * 5   #正电压端电位信号       V_+  乘以5以恢复量程
        self.V_f            = self.data[:, 9] * 5   #悬浮电位信号           V_f  乘以5以恢复量程
        self.Te             = self.electron_temperature()
        self.ne             = self.electron_density()

    def electron_temperature(self):
        e  = 1.602e-19
        k  = 1.38e-23
        Te = e * (self.V_plus - self.V_f) / (k * np.log(2))
        return Te

    def electron_density(self):
        C  = 0.909e19
        ne = self.V_ion * C
        return ne
    
    def magnetic_field_integral(self,B_data_1d,delta_time,s_eff,B0 = 0):
        """
        对磁场信号积分
        考虑物理意义，磁重联开始与结束后，磁场的诊断值应该为0（不考虑等离子体本身产生的背景磁场）
        """
        B_integral = np.array([])
        B_integral = np.append(B_integral,B0-delta_time*B_data_1d[0]/s_eff)
        B_int      = B0
        for i in range(1,len(B_data_1d)):
            #B_int = -np.sum(B_data_1d[:i]*delta_time/s_eff)
            B_int += -(B_data_1d[i-1]+B_data_1d[i])*delta_time/(2*s_eff)
            B_integral = np.append(B_integral,B_int)
        return B_integral
    
    def magnetic_signal_line_average(self,B_data_1,B_data_2):
        """
        平均双股线上的反向信号
        """
        return np.mean([B_data_1,-B_data_2],axis=0)

class DataProcess:
    """
    数据处理类，对指定一维数据创建时间窗，对齐时间测量误差
    """
    def __init__(self,
                 data_origin,
                 find_kneed_error = 0.006,
                 window_step = [100,500]
                 ):
        self.s = find_kneed_error
        self.before_kneed = window_step[0]
        self.after_kneed = window_step[1]

        self.kneed_point = self.find_kneed_arg(data_origin)
        self.kneed_d_point  = len(data_origin)-self.find_kneed_arg(data_origin[::-1])

        self.window_index = [self.kneed_point - self.before_kneed,
                             self.kneed_point + self.after_kneed]
        self.data_window = self.snap_data(data_origin, self.kneed_point,
                                          self.before_kneed, self.after_kneed)

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
    
    def find_kneed_down_args(self,data1d):
        """
        找到下降拐点的位置
        data1d: 一维数据
        s: 为阈值
        return 拐点的位置索引
        """
        diff = np.diff(data1d)
        point = np.where(diff > 0.004)
        point = point[0][0]
        return point

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

    def position_average_datas(self,datas):
        """
        对N维位置数据进行平均，采用np.mean的形式，可拓展
        datas: 高维位置数据
        """
        data_average = np.mean(datas)
        return data_average

def find_pathnum(data_path = './data/y3x3/'):
    '''得到mkr_path中y与x间的数字'''
    indy = data_path.find('y')
    indx = data_path.find('x')
    indz = data_path.find('/',indx)
    return data_path[indy:indz]


def position_average(mkr_path = './data/y4x15/'):
    Te = np.array([])
    ne = np.array([])
    Bxmi = np.array([])
    Bymi = np.array([])
    I_driven = np.array([])
    for i in np.arange(1,11,
                       #desc = 'Processing data under:'+mkr_path
                       ):
        data_path = mkr_path + 'N' +str(i) + '.mat'
        MRD = Magnetic_Reconnection_Diagnostic(data_path,
                                               B_Diagnostic=True)
        '''依据驱动电流找到磁重联起始点，创建时间窗'''
        DP = DataProcess(MRD.I_drive,window_step=[100,700])
        I_driven_window = DP.data_window
        window_index    = DP.window_index 

        '''对其它参数创建时间窗'''
        '静电探针获取计算的参数'
        Te_window = MRD.Te[window_index[0]:window_index[1]]
        ne_window = MRD.ne[window_index[0]:window_index[1]]
       

        '''磁探针信号需要重新处理，这里，认为再驱动电流启动前后，系统中的场忽略不计'''
        Bx = MRD.B_x_mean
        By = MRD.B_y_mean
        start_point = DP.kneed_point
        end_point = DP.kneed_d_point
        Bx[0:start_point] = 0
        Bx[end_point:]    = 0
        By[0:start_point] = 0
        By[end_point:]    = 0
        B_x_m_ig = MRD.magnetic_field_integral(Bx,MRD.delta_t,s_eff=3.001691131e-4)
        B_y_m_ig = MRD.magnetic_field_integral(By,MRD.delta_t,s_eff=4.550113814e-4)

        '磁探针获取计算的参数'
        Bxmi_window = B_x_m_ig[window_index[0]:window_index[1]]
        Bymi_window = B_y_m_ig[window_index[0]:window_index[1]]

        '''加入数组'''
        I_driven= np.append(I_driven, I_driven_window).reshape(-1, len(I_driven_window))
        Te      = np.append(Te, Te_window).reshape(-1, len(Te_window))
        ne      = np.append(ne, ne_window).reshape(-1, len(ne_window))
        Bxmi    = np.append(Bxmi,Bxmi_window).reshape(-1, len(Bxmi_window))
        Bymi    = np.append(Bymi,Bymi_window).reshape(-1, len(Bymi_window))


    '''对每个点的数据进行平均'''
    Te_average = np.mean(Te, axis=0)
    ne_average = np.mean(ne, axis=0)
    Bxmi_average = np.mean(Bxmi,axis=0)
    Bymi_average = np.mean(Bymi,axis=0)

    I_driven_average = np.mean(I_driven, axis=0)
    
    time = np.arange(0, len(I_driven_average), 1)



    scio.savemat('./dataall/data_processed/'+ find_pathnum(mkr_path) +'.mat',
                                               {'Te_average':Te_average,
                                                'ne_average':ne_average,
                                                'I_driven_average':I_driven_average,
                                                'Bxmi_average':Bxmi_average,
                                                'Bymi_average':Bymi_average,
                                                #'time':time
                                                })
    
    '''
    plt.figure(figsize=(20,5))
    for i in range(0,10):
        #plt.plot(time,Te[i])
        plt.plot(time,I_driven[i])
    #plt.plot(time, Te_average*1e13, label='Te')
    #plt.plot(time, ne_average, label='ne')
    #plt.plot(time, I_driven_average*2e17, label='I_driven')
    #plt.plot(time,Bxmi_average*2e20,label = 'Bxmi')
    #plt.plot(time,Bymi_average*2e20,label = 'Bymi')
    plt.legend()
    plt.show()'''
    
    
    
#position_average(mkr_path = './data/y4x15/N')

def per_processing_all_data(mutiprocess = False):
    if mutiprocess == True:
        import multiprocessing as mp
        from multiprocessing import RLock
        '''
        生成包含d1,d2的path的列表
        '''
        path_list = []
        for d1 in np.arange(1,14):
            for d2 in np.arange(1,26):
                data_path = './dataall/dataAll/y'+str(d1)+'x'+str(d2)+'/'
                path_list.append(data_path)
        '''
        多进程处理
        '''
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


if __name__ == '__main__':
    st = time.time()
    per_processing_all_data(mutiprocess = True)
    stt = time.time()
    print(stt-st)
    #position_average()