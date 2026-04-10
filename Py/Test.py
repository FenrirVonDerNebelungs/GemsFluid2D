import CTrans as ct
import numpy as np
import matplotlib.pyplot as plt

g_grid_width=4*16
g_blow_factor = 24
g_num_grid_standalone_frames = 11
g_num_jacobi_frames = 10
g_num_expanded_frames = 2
g_delta_x=1e-3
g_delta_t=1e-3

class drawFrame:
    def __init__(self, frame_index, filename, num_headers_per_frame, framesize):
        self.grid_wh = g_grid_width
        self.grid_exp_wh = g_grid_width*g_blow_factor
        self.frame_index=frame_index

    def updateFrameIndex(self, frame_index):
        self.frame_index = frame_index

    def convertTo2D(self, raw_array, exp_factor=1):
        width = g_grid_width*exp_factor
        height = width
        V_array = np.zeros((width, height))
        for j in range(height):
            for i in range(width):
                V_array[i,j]=raw_array[j*width + i]
        return V_array


    def draw_U_vector(self, Ux, Uy, graph_title):
        x,y = np.meshgrid(np.arrange(g_delta_x, g_grid_width*g_delta_x, g_delta_x), np.arrange(g_delta_x, g_grid_width*g_delta_x, g_delta_x))
        Ux_2D = self.convertTo2D(Ux)
        Uy_2D = self.convertTo2D(Uy)
        Mag = np.sqrt(Ux_2D**2 + Uy_2D**2)
        plt.figure(figsize=(7,7))
        plt.quiver(x,y,Ux_2D, Uy_2D, Mag, cmap='viridis')
        plt.colorbar(label='Magnitude')

        plt.title(graph_title)
        plt.grid(True)
        plt.axis('equal')

        plt.show()

    def draw_scalar(self, p, graph_title, exp_factor=1):
        delta=g_delta_x/(float(exp_factor))
        mesh_start = delta/2.0 #inclusive, defaults to zero if not specified
        mesh_stop = g_grid_width*g_delta_x + delta/2.0 #exclusive not in array
        mesh_step = delta
        x,y = np.meshgrid(np.arrange(mesh_start, mesh_stop, mesh_step), np.arrange(mesh_start, mesh_stop, mesh_step))
        p_2D = self.convertTo2D(p, exp_factor)
        fig=plt.figure(figsize=(7,7))
        ax=plt.axes(projection='3d')
        surf = ax.plot_surface(x,y,p_2D,cmap='coolwarm')
        ax.set_title(graph_title)
        plt.show()

    def draw_scalar_and_points(self, U_exp, U, graph_title, exp_factor):
        fig=plt.figure(figsize=(7,7))
        ax = fig.add_subplot(111, projection='3d')

        delta=g_delta_x/(float(exp_factor))
        mesh_start = delta/2.0 #inclusive, defaults to zero if not specified
        mesh_stop = g_grid_width*g_delta_x + delta/2.0 #exclusive not in array
        mesh_step = delta
        x_exp,y_exp = np.meshgrid(np.arrange(mesh_start, mesh_stop, mesh_step), np.arrange(mesh_start, mesh_stop, mesh_step))
        U_exp_2D = self.convertTo2D(U_exp, exp_factor)

        surf = ax.plot_surface(x_exp,y_exp,U_exp_2D,cmap='coolwarm')

        mesh_start = g_delta_x/2.0
        mesh_stop = g_grid_width*g_delta_x + g_delta_x/2.0
        mesh_step = g_delta_x
        x,y = np.meshgrid(np.arrange(mesh_start, mesh_stop, mesh_step), np.arrange(mesh_start, mesh_stop, mesh_step))
        U_2D = self.convertTo2D(U)
        ax.scatter(x, y, U_2D, color='green', s=50)

        ax.set_title(graph_title)
        plt.show()

    def draw(self, header_stack, data_stack):
        for i in range(len(header_stack)):
            header = header_stack[i]
            label = ct.header2[ct.getDataLabel(header)]
            axis = ct.header3[ct.getDataAxis(header)]
            start_end_label = ct.header4[ct.getDataStartEndCode(header)]
            match label:
                case 'U':
                    expansion_factor = ct.getExpansionFactor(header)
                    if expansion_factor==1 and axis=='X':
                        Ux = data_stack[i]
                        Uy = data_stack[i+1]
                        graph_title = "U at Frame: "+str(self.frame_index) + " step:  " +start_end_label
                        self.draw_U_vector(Ux, Uy, graph_title)
                    else:
                        U_exp = data_stack[i]
                        graph_title = "U "+axis +"  expanded"
                        self.draw_scalar(U_exp, graph_title, expansion_factor)
                case'relPos':
                    expansion_factor = ct.getExpansionFactor(header)
                    U_exp = data_stack[i-2]
                    U_back = data_stack[i]
                    graph_title = "Back traced U: "+axis
                    self.draw_scalar_and_points(U_exp, U_back, graph_title, expansion_factor)
                case 'W':
                    expansion_factor = ct.getExpansionFactor(header)
                    if expansion_factor==1:
                        Wx = data_stack[i]
                        Wy = data_stack[i+1]
                        graph_title = "W at Frame: "+str(self.frame_index) 
                        self.draw_U_vector(Wx, Wy, graph_title)
                case 'DivW':
                    DivW = data_stack[i]
                    self.draw_scalar(DivW, 'Div W')
                case 'jacobi_frame':
                    frame_i = ct.getJacobiFrame(header)
                    pressure = data_stack[i]
                    frame_title = 'P jacobi frame: '+ str(frame_i)
                    self.draw_scalar(pressure, frame_title)
                case 'P':
                    P = data_stack[i]
                    self.draw_scalar(pressure, label)
                case 'gradP':
                    x=1
                case 'DivU':
                    DivU = data_stack[i]
                    self.draw_scalar(DivU, label)
                case _:
                    xxx='unmatched'
            
        return 0

class Test:

    def __init__(self, filename = '../Dat/frames.dat'):#'../Dat/test.dat'):#'../Dat/frames.dat'):
        grid_stream_len = g_grid_width*g_grid_width
        grid_exp_stream_len = grid_stream_len*g_blow_factor*g_blow_factor
        total_grid_stream_len = grid_stream_len*(11+g_num_jacobi_frames)+grid_exp_stream_len*2
        total_num_stream_headers = (11+g_num_jacobi_frames)+2
        header_len_bytes = ct.CTrans.header_len
        total_grid_stream_len_bytes = total_grid_stream_len*8
        total_headers_len_bytes = header_len_bytes*total_num_stream_headers
        self.framesize = total_grid_stream_len_bytes + total_headers_len_bytes
        self.filename = filename
        self.num_headers = total_num_stream_headers
        self.C_Trans = ct.CTrans(filename, total_num_stream_headers,self.framesize)
        self.Draw = drawFrame(0,filename,total_num_stream_headers,self.framesize)

    def testRun(self):
        dat, file_offset = self.C_Trans.readDdat(4,0)
        print(dat)

    def Run(self):
        len_returned = self.num_headers
        frame_cnt=0
        file_offset=0
        while len_returned>1:
            header_stack, data_stack, new_file_offset = self.C_Trans.readFrame(file_offset)
            file_offset=new_file_offset
            len_returned = len(header_stack)
            self.Draw.updateFrameIndex(frame_cnt)
            frame_cnt += 1
            if(len_returned==self.num_headers):
                self.Draw.draw(header_stack, data_stack)
            
testInst = Test()
testInst.Run()
print("finished")