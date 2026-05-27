import CTrans as ct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def findMaxMin(raw_array):
    max_val=0.0
    min_val=0.0
    max_val = np.max(raw_array)
    min_val = np.min(raw_array)
    return max_val, min_val

def convertDoubleToFloatArray(double_array):
    np_array = np.array(double_array)
    return np_array.astype(np.float32)

def convertTo2D(raw_array, width, height):
    V_array = raw_array.reshape(height, width) 
    return V_array

def convertDoubleTo2DFloat(double_array, width, height):
    float_array = convertDoubleToFloatArray(double_array)
    return convertTo2D(float_array,width,height)

class ViewFrames:
    def __init__(self, filename='../Dat/frames.dat', delta_x=1.0e-3):
        self.delta_x = delta_x
        self.filename=filename
        self.max_frames=1000
        self.frame_dt = 1000
        self.C_Trans = ct.CTrans(filename)
        # dim of P,Ux,Uy,Dy frames * 2D array of width x height
        self.P = []
        self.Ux = []
        self.Uy = []
        self.Dye = []
        self.width = -1
        self.height = -1
        self.error_flag = False
    
        figP, axP =  plt.subplots(figsize=(6,6))
        self.figP = figP
        self.axP = axP
        self.imP = None

        figDye,axDye = plt.subplots(figsize=(6,6))
        self.figDye = figDye
        self.axDye = axDye
        self.imDye=None
        figU,axU = plt.subplots(figsize=(6,6))
        self.figU=figU
        self.axU =axU
        self.imU = None
        #dim of x, y 2D array of width x height
        self.x =[]
        self.y =[]
        self.x_sparce=[]
        self.y_sparce=[]
        self.step=1
        
        self.aniP=None
        self.aniDye=None
        self.aniU=None
    
    def setFrame_dt(self,dt):
        if dt>=1 and dt<=10:
            self.frame_dt = 100*dt
    
    def runAni(self):
        self.assembleFrames()
        self.firstFrame()
        if not self.error_flag:
            num_frames = len(self.P)
            self.aniP = FuncAnimation(self.figP, self.update_P, frames=num_frames, interval=self.frame_dt, blit=True)
            self.figP.show()
            num_frames=len(self.Dye)
            self.aniDye = FuncAnimation(self.figDye, self.update_Dye, frames=num_frames, interval=self.frame_dt, blit=True)
            self.figDye.show()
            num_frames = len(self.Ux)
            self.aniU = FuncAnimation(self.figU, self.update_U,frames=num_frames, interval=self.frame_dt, blit=True)
            self.figU.show()
    
    def firstFrame(self):
        if(self.error_flag):
            return
        gwidth=self.width
        gheight = self.height
        mesh_start = self.delta_x
        mesh_stop_x = (gwidth+1)*self.delta_x
        mesh_stop_y = (gheight+1)*self.delta_x
        mesh_step = self.delta_x
        x_quiv,y_quiv = np.meshgrid(np.arange(mesh_start, mesh_stop_x, mesh_step), np.arange(mesh_start, mesh_stop_y, mesh_step))
        self.step=1
        if(gwidth>64):
            self.step=int(4*gwidth/128)
        self.x_sparse = x_quiv[::self.step, ::self.step]
        self.y_sparse = y_quiv[::self.step, ::self.step]
        Ux_sparse = self.Ux[0][::self.step, ::self.step]
        Uy_sparse = self.Uy[0][::self.step, ::self.step]
        Mag = np.sqrt(Ux_sparse**2 + Uy_sparse**2)
        self.imU = self.axU.quiver(self.x_sparse,self.y_sparse,Ux_sparse, Uy_sparse, Mag, angles='xy', scale_units='xy', scale=50.0, cmap='viridis', animated=True)
    
        mesh_start = self.delta_x/2.0 #inclusive, defaults to zero if not specified
        mesh_stop = gwidth*self.delta_x + self.delta_x/2.0 #exclusive not in array
        self.x,self.y = np.meshgrid(np.arange(mesh_start, mesh_stop, mesh_step), np.arange(mesh_start, mesh_stop, mesh_step))
        self.imP = self.axP.imshow(self.P[0], cmap='jet', interpolation='gaussian', origin='lower', aspect='equal', animated=True)
        self.imDye = self.axDye.imshow(self.Dye[0], cmap='twilight', interpolation='gaussian', origin='lower', aspect='equal', animated=True)
        self.imP.set_extent([self.x.min(), self.x.max(), self.y.min(), self.y.max()])
        self.imDye.set_extent([self.x.min(), self.x.max(), self.y.min(), self.y.max()])

    def assembleFrames(self):
        len_returned = 2
        num_frames_read=0
        print("\n \n Reading file: ",end="")
        while len_returned>1 and num_frames_read<=self.max_frames:
            header_stack, data_stack, new_file_offset = self.C_Trans.readFrame()
            num_frames_read+=1
            len_returned =len(header_stack)
            if(len_returned>=1):
                for i_stack in range(len(header_stack)):
                    header = header_stack[i_stack]
                    grid_width = ct.getGridWidth(header)
                    grid_height = ct.getGridHeight(header)
                    if grid_width!=self.width or grid_height!=self.height:
                        if self.width > 0 or self.height >0:
                            print("Error reading grid dimensions")
                            self.error_flag=True
                            break
                        else:
                            self.width=grid_width 
                            self.height = grid_height
                    label = ct.header2[ct.getDataLabel(header)]
                    axis = ct.header3[ct.getDataAxis(header)]
                    match label:
                        case 'U':
                            if axis == 'X':
                                red_datUx = convertDoubleTo2DFloat(data_stack[i_stack], self.width,self.height)
                                self.Ux.append(red_datUx)
                            elif axis == 'Y':
                                red_datUy = convertDoubleTo2DFloat(data_stack[i_stack],self.width,self.height)
                                self.Uy.append(red_datUy)
                        case 'P':
                            red_dat = convertDoubleTo2DFloat(data_stack[i_stack],self.width,self.height)
                            self.P.append(red_dat)
                        case 'Dye':
                            red_dat = convertDoubleTo2DFloat(data_stack[i_stack],self.width,self.height)
                            self.Dye.append(red_dat)
            if(self.error_flag):
                break
            if num_frames_read%10==0:
                print('*',end="")
        print('!')
        return

    def update_P(self, frame):
        #frame = i_stack
        P_dat = self.P[frame]
        max,min = findMaxMin(P_dat)
        self.imP.set_data(P_dat) 
        self.imP.set_clim(vmin=min, vmax=max)
        return [self.imP]

    def update_Dye(self,frame):
        Dye_dat = self.Dye[frame]
        max,min=findMaxMin(Dye_dat)
        self.imDye.set_data(Dye_dat)
        self.imDye.set_clim(vmin=min, vmax=max)
        return [self.imDye]

    def update_U(self,frame):
        Ux_dat = self.Ux[frame]
        Uy_dat = self.Uy[frame]
        Ux_sparce = Ux_dat[::self.step, ::self.step]
        Uy_sparce = Uy_dat[::self.step, ::self.step]
        Mag = np.sqrt(Ux_sparce**2 + Uy_sparce**2)
        self.imU.set_UVC(Ux_sparce, Uy_sparce,Mag)
        return [self.imU]

ViewFrames_inst = ViewFrames()
dt_str = input("Set speed: 10 slow,  1 fast: ")
dt = int(dt_str)
ViewFrames_inst.setFrame_dt(dt)
ViewFrames_inst.runAni()
exit_val='n'
while exit_val!='y':
    exit_val = input("y to exit: ")
    print(exit_val)
print("finished")