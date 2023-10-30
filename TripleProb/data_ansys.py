import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import scipy.io as scio
from matplotlib import gridspec as gp

import multiprocessing as mp
from tqdm import tqdm

#来自自建文件
from data_process import DataProcess

#指定ffmpeg
plt.rcParams['animation.ffmpeg_path']='C:\\Users\\USTC\\Documents\\ffmpeg-6.0-essentials\\bin\\ffmpeg.exe'

class DrawAndAnimation():
    def __init__(self,
        
                 ):
        self.dt = 1/2e6
        self.Tensor_generate()
    
    def Tensor_generate(self):
        '''构建每个数据的Tensor'''

        '''为了获取数据维度'''
        data_path = './dataall/data_processed/y6x12.mat'
        data      = scio.loadmat(data_path)
        timett    = len(data['Te_average'][0])

        I  = data['I_driven_average'][0]
        DP          = DataProcess(I)
        start_point = DP.before_kneed
        end_point   = DP.kneed_d_point

        Bx_tensor = np.zeros((timett,13,25))
        By_tensor = np.zeros((timett,13,25))
        Bt_tensor = np.zeros((timett,13,25))

        Te_tensor = np.zeros((timett,13,25))
        ne_tensor = np.zeros((timett,13,25))

        for d1 in np.arange(1,14):
            for d2 in np.arange(1,26):
                data_path = './dataall/data_processed/y'+str(d1)+'x'+str(d2)+'.mat'
                data = scio.loadmat(data_path)
                
                Bx = data['Bxmi_average'][0]
                #Bx[0:start_point] = 0
                Bx[end_point:] = 0
                Bx_tensor[:,d1-1,d2-1] = Bx

                By = data['Bymi_average'][0]
                #By[0:start_point] = 0
                By[end_point:] = 0
                By_tensor[:,d1-1,d2-1] = By

                Bt = np.sqrt(Bx**2+By**2)
                Bt_tensor[:,d1-1,d2-1] = Bt

                Te = data['Te_average'][0]
                Te_tensor[:,d1-1,d2-1] = Te

                ne = data['ne_average'][0]
                ne_tensor[:,d1-1,d2-1] = ne
        
        #将Te在x方向反转，为什么要反转呢？是否由于streamplot的原因？
        Te_tensor = Te_tensor[:,:,::-1]
        ne_tensor = ne_tensor[:,:,::-1]
        #return ne_tensor,Te_tensor,Bx_tensor,By_tensor,Bt_tensor 
        self.ne_tensor = ne_tensor
        self.Te_tensor = Te_tensor
        self.Bx_tensor = Bx_tensor
        self.By_tensor = By_tensor
        self.Bt_tensor = Bt_tensor
        self.start_point = start_point
        self.end_point = end_point
        self.I         = I

    def Animation_framework(self,
                            compont,
                            compont_name,
                            before_start_time = 50,
                            after_end_time = 50,
                            BFL=True,
                            save=True):
        cmax = np.max(compont)
        cmin = np.min(compont)

        x1d = np.arange(0,25)
        x2d = np.arange(0,13)
        x1,y1 = np.meshgrid(x1d,x2d)

        "生成streanline的startpoint数组"
        start_points = []
        for j in np.linspace(1,23,50):
            start_points.append([j,0])
            start_points.append([j,12])

        fig,ax = plt.subplots(figsize=(12,5))
        cf = plt.contourf(compont[0,:,:],150,cmap = 'plasma',levels=np.linspace(cmin,cmax,150))
        fig.colorbar(cf)

        def update(frame):
            frame = int(frame+self.start_point-before_start_time)
            ax.clear()
            ax.contourf(compont[frame,:,:],150,levels=np.linspace(cmin,cmax,150),cmap='plasma')
            if BFL == True:
                ax.streamplot(x=x1,y=y1,u=self.Bx_tensor[frame,y1,x1],v=-self.By_tensor[frame,y1,x1],
                    color=self.Bt_tensor[frame,y1,x1],
                    cmap='cool',
                    broken_streamlines=False,
                    linewidth=0.5,
                    arrowsize = 0.5,
                    start_points = start_points,
                    )
            ax.set_title('{:}@t = {:.5f}ms'.format(compont_name,(-self.start_point+frame)*self.dt*1e3))
            if frame>=self.start_point and frame<=self.end_point:
                ax.set_title('{:}@t = {:.5f}ms:Reconnection Begun'.format(compont_name,(-self.start_point+frame)*self.dt*1e3))   

        frame_total = self.end_point-self.start_point+before_start_time+after_end_time
        ani = FuncAnimation(fig,
                    func = update,
                    frames = frame_total,
                    interval = 0.005
        )

        if save == True:
            '''创建进度条'''
            pos = np.random.randint(4) #取巧的做法，理论上有1/4的概率使得进度条重叠
            pbar = tqdm(total=frame_total,position=pos,desc="Saving {:}".format(compont_name))
            def progress_callback(current_frame,total_frame):
                pbar.update(1)
                
            print('\n{:} now Saving...\n'.format(compont_name))
            ani.save('{:} with BFL={:}.mp4'.format(compont_name,BFL),
                     fps = 60,
                     writer='ffmpeg',
                     dpi = 330,
                     progress_callback= progress_callback
                     )  
        else:
            plt.show()
        
    def BigAnimation(self,
                     ways = 'x1y4',
                     before_start_time = 50,
                     after_end_time = 50):
        Tmax = np.max(self.Te_tensor)
        Tmin = np.min(self.Te_tensor)

        nmax = np.max(self.ne_tensor)
        nmin = np.min(self.ne_tensor)

        Bmax = np.max(self.Bt_tensor)
        Bmin = np.min(self.Bt_tensor)

        x1d = np.arange(0,25)
        x2d = np.arange(0,13)
        x1,y1 = np.meshgrid(x1d,x2d)

        "生成streanline的startpoint数组"
        start_points = []
        for j in np.linspace(1,23,50):
            start_points.append([j,0])
            start_points.append([j,12])

        if ways == 'x3y2':
            fig,([ax,ax1],[ax2,ax3]) = plt.subplots(2,2,figsize=(24,10))
            fig = plt.figure(figsize=(28,6))
            plt.subplots_adjust(left=0.02, bottom=0.11, right=0.97, top=0.88,wspace= 0.05)
            gs = gp.GridSpec(2,3,height_ratios=[2,1],width_ratios=[1,1,1])
            ax = plt.subplot(gs[0,0])
            ax1 = plt.subplot(gs[0,1])
            ax2 = plt.subplot(gs[0,2])
            ax3 = plt.subplot(gs[1,:])
        elif ways == 'x2y2':
            fig,([ax,ax1],[ax2,ax3]) = plt.subplots(2,2,figsize=(24,10))
            fig = plt.figure(figsize=(28,6))
            plt.subplots_adjust(left=0.02, bottom=0.11, right=0.97, top=0.88,wspace= 0.05)
            gs = gp.GridSpec(2,2,height_ratios=[1,1],width_ratios=[1,1])
            ax = plt.subplot(gs[0,0])
            ax1 = plt.subplot(gs[0,1])
            ax2 = plt.subplot(gs[1,0])
            ax3 = plt.subplot(gs[1,1])
        elif ways == 'x1y4':
            fig = plt.figure(figsize=(12,12))
            plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.96,wspace= 0.05)
            gs = gp.GridSpec(4,1,height_ratios=[1,1,1,0.5])
            ax = plt.subplot(gs[0,:])
            ax1 = plt.subplot(gs[1,:])
            ax2 = plt.subplot(gs[2,:])
            ax3 = plt.subplot(gs[3,:])






        cf = plt.contourf(self.Te_tensor[0,:,:],150,cmap = 'plasma',levels=np.linspace(Tmin,Tmax,150))
        cf2 = plt.contourf(self.ne_tensor[0,:,:],150,cmap = 'plasma',levels=np.linspace(nmin,nmax,150))
        cf3 = plt.contourf(self.Bt_tensor[0,:,:],150,cmap = 'plasma',levels=np.linspace(Bmin,Bmax,150))

        frame_total = self.end_point-self.start_point+before_start_time+after_end_time
        
        #cf4=ax3.plot(int(self.start_point-before_start_time)+np.arange(0,frame_total),
         #              self.I[self.start_point-before_start_time:self.end_point+after_end_time])
        ax3.set_ylim(1.1*np.min(self.I),1.1*np.max(self.I))
        ax3.set_xlim(self.start_point-before_start_time,self.end_point+after_end_time)
        ax3.text(self.start_point+10,self.I[self.start_point+10],'')
        ax3.set_position([ax3.get_position().x0, ax3.get_position().y0, 0.8*ax3.get_position().width, ax3.get_position().height])         
        plt.colorbar(cf,ax=ax)
        plt.colorbar(cf2,ax=ax1)
        plt.colorbar(cf3,ax=ax2)
        

        def update(frame):
            frame = int(frame+self.start_point-before_start_time)
            ax.clear()
            ax1.clear()
            ax2.clear()
            ax3.clear()
            ax.contourf(self.Te_tensor[frame,:,:],150,levels=np.linspace(Tmin,Tmax,150),cmap='plasma')

            #绘制密度部分
            ax1.contourf(self.ne_tensor[frame,:,:],150,levels=np.linspace(nmin,nmax,150),cmap='plasma')

            #绘制磁场部分
            ax2.contourf(self.Bt_tensor[frame,:,:],150,levels=np.linspace(Bmin,Bmax,150),cmap='plasma')

            #绘制电流曲线曲线部分
            ax3.scatter(frame,self.I[frame],c='r')
            ax3.plot(int(self.start_point-before_start_time)+np.arange(0,frame_total),
                       self.I[self.start_point-before_start_time:self.end_point+after_end_time])
            ax3.axvline(frame,c='r',lw=1)
            ax3.axhline(self.I[frame],c='r',lw=1)
            ax3.text(self.start_point-before_start_time-10,self.I[frame],'I={:.5f}'.format(self.I[frame]))
            
            
            
            ax.streamplot(x=x1,y=y1,u=self.Bx_tensor[frame,y1,x1],v=-self.By_tensor[frame,y1,x1],
                color=self.Bt_tensor[frame,y1,x1],
                cmap='cool',
                broken_streamlines=False,
                linewidth=0.5,
                arrowsize = 0.5,
                start_points = start_points,
                )
            ax1.streamplot(x=x1,y=y1,u=self.Bx_tensor[frame,y1,x1],v=-self.By_tensor[frame,y1,x1],
                color=self.Bt_tensor[frame,y1,x1],
                cmap='cool',
                broken_streamlines=False,
                linewidth=0.5,
                arrowsize = 0.5,
                start_points = start_points,
                )
            ax2.streamplot(x=x1,y=y1,u=self.Bx_tensor[frame,y1,x1],v=-self.By_tensor[frame,y1,x1],
                color=self.Bt_tensor[frame,y1,x1],
                cmap='cool',
                broken_streamlines=False,
                linewidth=0.5,
                arrowsize = 0.5,
                start_points = start_points,
                )
            
            ax.set_title('{:}@t = {:.5f}ms'.format('Te/eV',(-self.start_point+frame)*self.dt*1e3))
            ax1.set_title('{:}@t = {:.5f}ms'.format('ne',(-self.start_point+frame)*self.dt*1e3))
            ax2.set_title('{:}@t = {:.5f}ms'.format('|B|',(-self.start_point+frame)*self.dt*1e3))
            ax3.set_title('{:}@t = {:.5f}ms'.format('I_driven/kA',(-self.start_point+frame)*self.dt*1e3))
            if frame>=self.start_point and frame<=self.end_point:
                ax.set_title('{:}@t = {:.5f}ms:Reconnection Begun'.format('Te/eV',(-self.start_point+frame)*self.dt*1e3)) 
                ax1.set_title('{:}@t = {:.5f}ms:Reconnection Begun'.format('ne',(-self.start_point+frame)*self.dt*1e3))
                ax2.set_title('{:}@t = {:.5f}ms:Reconnection Begun'.format('|B|',(-self.start_point+frame)*self.dt*1e3))
                ax3.set_title('{:}@t = {:.5f}ms:Reconnection Begun'.format('I_driven/kA',(-self.start_point+frame)*self.dt*1e3))  

        
        ani = FuncAnimation(fig,
                    func = update,
                    frames = frame_total,
                    interval = 0.005
        )
        plt.show()
        
        '''创建进度条'''
        pos = np.random.randint(4) #取巧的做法，理论上有1/4的概率使得进度条重叠
        pbar = tqdm(total=frame_total,position=pos,desc="Saving {:}".format('Te'))
        def progress_callback(current_frame,total_frame):
            pbar.update(1)
            
        print('\n{:} now Saving...\n'.format('Te'))
        ani.save('Bigwith41infps60.mp4',
                    fps = 60,
                    writer='ffmpeg',
                    dpi = 150,
                    progress_callback= progress_callback,
                    )  
        

def DrawTwoatOnce():
    '''1 实例化'''
    DA = DrawAndAnimation()
    Te = DA.Te_tensor
    ne = DA.ne_tensor
    Bx = DA.Bx_tensor
    By = DA.By_tensor
    Bt = DA.Bt_tensor

    pool = mp.Pool(4)
    pool.apply_async(DA.Animation_framework,
                    args=(Te,'Te',)
                    )
    pool.apply_async(DA.Animation_framework,
                    args=(ne,'ne',)
                    )
    pool.close()
    pool.join()
    """ 
   p1 = mp.Process(target=DA.Animation_framework,
                    args=(Te,'Te',)
                    )
    
    p2 = mp.Process(target=DA.Animation_framework,
                    args=(ne,'ne',)
                    )

     p1.start()
    p2.start()
    p1.join()
    p2.join()
    """



   

if __name__ == "__main__":
   DA = DrawAndAnimation()
   DA.BigAnimation()