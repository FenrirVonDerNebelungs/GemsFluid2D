import CTrans as ct
import numpy as np
import matplotlib.pyplot as plt
import TestCheck as tc


g_delta_x=1e-3
g_delta_t=1e-3

class drawFrame:
    def __init__(self, frame_index):
        self.grid_wh = [0,0]
        self.frame_index=frame_index
        self.jacobi_max=0.0
        self.jacobi_err_max =0.0

    def updateFrameIndex(self, frame_index):
        self.frame_index = frame_index
        
    def setGridWH(self, width, height):
        self.grid_wh = [width, height]

    def convertTo2D(self, raw_array, exp_factor=1):
        width = self.grid_wh[0]
        height = self.grid_wh[1]
        V_array = np.array(raw_array).reshape(height, width) #np.zeros((width, height))
        #for j in range(height):
        #    for i in range(width):
        #        V_array[i,j]=raw_array[j*width + i]
        return V_array

    def findMaxMin(self, raw_array):
        max_val=0.0
        min_val=0.0
        max_val = max(raw_array)
        min_val = min(raw_array)
        return max_val, min_val
    
    def draw_U_vector(self, Ux, Uy, graph_title):
        gwidth = self.grid_wh[0]
        gheight = self.grid_wh[1]
        x,y = np.meshgrid(np.arange(g_delta_x, (gwidth+1)*g_delta_x, g_delta_x), np.arange(g_delta_x, (gheight+1)*g_delta_x, g_delta_x))
        Ux_2D = self.convertTo2D(Ux)
        Uy_2D = self.convertTo2D(Uy)
        Mag = np.sqrt(Ux_2D**2 + Uy_2D**2)
        plt.figure(figsize=(12,12))
        plt.quiver(x,y,Ux_2D, Uy_2D, Mag, angles='xy', scale_units='xy', scale=1, cmap='viridis')
        plt.colorbar(label='Magnitude')

        plt.title(graph_title)
        plt.axis('equal')

        plt.show()

    def draw_scalar(self, p, graph_title, z_axis=-1, exp_factor=1):
        gwidth = self.grid_wh[0]
        delta=g_delta_x/(float(exp_factor))
        mesh_start = delta/2.0 #inclusive, defaults to zero if not specified
        mesh_stop = gwidth*g_delta_x + delta/2.0 #exclusive not in array
        mesh_step = delta
        x,y = np.meshgrid(np.arange(mesh_start, mesh_stop, mesh_step), np.arange(mesh_start, mesh_stop, mesh_step))
        p_2D = self.convertTo2D(p, exp_factor)
        fig=plt.figure(figsize=(7,7))
        ax=plt.axes(projection='3d')
        if z_axis > 0:
            ax.set_zlim(-z_axis, z_axis)
        surf = ax.plot_surface(x,y,p_2D,cmap='coolwarm')
        ax.set_title(graph_title)
        plt.show()

    def draw_scalar_and_points(self, U_exp, U, graph_title, exp_factor):
        gwidth = self.grid_wh[0]
        
        fig=plt.figure(figsize=(7,7))
        ax = fig.add_subplot(111, projection='3d')

        delta=g_delta_x/(float(exp_factor))
        mesh_start = delta/2.0 #inclusive, defaults to zero if not specified
        mesh_stop = gwidth*g_delta_x + delta/2.0 #exclusive not in array
        mesh_step = delta
        x_exp,y_exp = np.meshgrid(np.arange(mesh_start, mesh_stop, mesh_step), np.arange(mesh_start, mesh_stop, mesh_step))
        U_exp_2D = self.convertTo2D(U_exp, exp_factor)

        surf = ax.plot_surface(x_exp,y_exp,U_exp_2D,cmap='coolwarm')

        mesh_start = g_delta_x/2.0
        #mesh_stop = self.grid_wh*g_delta_x + g_delta_x/2.0
        mesh_step = g_delta_x
        x,y = np.meshgrid(np.arange(mesh_start, mesh_stop, mesh_step), np.arange(mesh_start, mesh_stop, mesh_step))
        U_2D = self.convertTo2D(U)
        ax.scatter(x, y, U_2D, color='green', s=50)

        ax.set_title(graph_title)
        plt.show()

    def drawColor2D(self, p, graph_title, max_val, min_val):
        gwidth = self.grid_wh[0]
        delta=g_delta_x
        mesh_start = delta/2.0 #inclusive, defaults to zero if not specified
        mesh_stop = gwidth*g_delta_x + delta/2.0 #exclusive not in array
        mesh_step = delta
        x,y = np.meshgrid(np.arange(mesh_start, mesh_stop, mesh_step), np.arange(mesh_start, mesh_stop, mesh_step))
        p_2D = self.convertTo2D(p)
        fig=plt.figure(figsize=(7,7))
        plt.imshow(p_2D, 
                   cmap='jet', 
                   interpolation='gaussian', 
                   origin='lower', 
                   aspect='equal',
                   vmax=max_val,
                   vmin=min_val, 
                   extent=[x.min(), x.max(), y.min(), y.max()])
        plt.title(graph_title)
        plt.colorbar()
        plt.show()

    def draw(self, header_stack, data_stack):
        for i in range(len(header_stack)):
            header = header_stack[i]
            grid_width = ct.getGridWidth(header)
            grid_height = ct.getGridHeight(header)
            self.setGridWH(grid_width, grid_height)
            label = ct.header2[ct.getDataLabel(header)]
            axis = ct.header3[ct.getDataAxis(header)]
            start_end_label = ct.header4[ct.getDataStartEndCode(header)]
            jacobi_frame_label='jacobi'
            jacobi_fill_label = 'fill'
            if (not jacobi_frame_label in start_end_label) and (not jacobi_fill_label in start_end_label):
                match label:
                    case 'U':
                        if axis=='X':
                            Ux = data_stack[i]
                            Uy = data_stack[i+1]
                            graph_title = "U at Frame: "+str(self.frame_index) + " step:  " +start_end_label
                            self.draw_U_vector(Ux, Uy, graph_title)
                    case'relPos':
                        if axis=='X':
                            rel_x = data_stack[i]
                            rel_y = data_stack[i+1]
                            graph_title = "Back traced U after delta_t: "
                            self.draw_U_vector(rel_x, rel_y, graph_title)
                    case 'W':
                        if axis=='X':
                            Wx = data_stack[i]
                            Wy = data_stack[i+1]
                            graph_title = "W at Frame: "+str(self.frame_index)+ " step: "+start_end_label 
                            self.draw_U_vector(Wx, Wy, graph_title)
                    case 'DivW':
                        DivW = data_stack[i]
                        self.draw_scalar(DivW, 'Div W')
                    case 'P':
                        pressure = data_stack[i]
                        graph_title = label + " at Frame: "+str(self.frame_index) + " step: "+start_end_label
                        self.draw_scalar(pressure, label)
                    case 'DivU':
                        DivU = data_stack[i]
                        self.draw_scalar(DivU, label)
                    case _:
                        xxx='unmatched'
            
        return 0
    
    def get_W_DivW_P_GradP_limits(self, header_stack, data_stack):
        W_min = 0.0
        W_max = 0.0
        DivW_min = 0.0
        DivW_max = 0.0
        P_min = 0.0
        P_max = 0.0
        gradP_min=0.0
        gradP_max=0.0
        
        for i in range(len(header_stack)):
            header = header_stack[i]
            label = ct.header2[ct.getDataLabel(header)]
            start_end_label = ct.header4[ct.getDataStartEndCode(header)]
            jacobi_frame_label = 'jacobi'
            if jacobi_frame_label in start_end_label:
                scalar_val = data_stack[i]
                cur_max, cur_min = self.findMaxMin(scalar_val)
                if label == 'W':
                    if(cur_max>W_max):
                        W_max = cur_max
                    if(cur_min<W_min):
                        W_min = cur_min
                elif label == 'DivW':
                    if(cur_max>DivW_max):
                        DivW_max = cur_max
                    if(cur_min<DivW_min):
                        DivW_min = cur_min
                elif label == 'P':
                    if(cur_max>P_max):
                        P_max=cur_max
                    if(cur_min<P_min):
                        P_min=cur_min
            else:
                if label=='gradP':
                    if(cur_max>gradP_max):
                        gradP_max = cur_max
                    if(cur_min<gradP_min):
                        gradP_min = cur_min
        return W_min,W_max,DivW_min,DivW_max,P_min,P_max,gradP_min,gradP_max
    
    def drawJacobi(self, header_stack, data_stack):
        W_min, W_max, DivW_min, DivW_max, P_min, P_max, gradP_min, gradP_max = self.get_W_DivW_P_GradP_limits(header_stack,data_stack)                                        
        for i in range(len(header_stack)):
            header = header_stack[i]
            grid_width = ct.getGridWidth(header)
            grid_height = ct.getGridHeight(header)
            self.setGridWH(grid_width, grid_height)
            numLoops = ct.getLoopsRun(header)
            label = ct.header2[ct.getDataLabel(header)]
            axis_label = ct.header3[ct.getDataAxis(header)]
            start_end_label = ct.header4[ct.getDataStartEndCode(header)]
            jacobi_frame_label = 'jacobi'
            if jacobi_frame_label in start_end_label:
                frame_i = ct.getJacobiFrame(header)
                stack_i = ct.getJacobiStackI(header)
                frame_title = label + "  stack: " +str(stack_i)
                scalar_val = data_stack[i]
                if start_end_label=='jacobi_stack_fill':
                    frame_title += ' fill'
                    cur_max=0.0
                    cur_min=0.0
                    if axis_label != 'Scalar':
                        frame_title += ' '+axis_label
                        cur_max = W_max
                        cur_min = W_min
                    elif label == 'DivW':
                        cur_max = DivW_max
                        cur_min = DivW_min
                    else:
                        cur_max=P_max
                        cur_min=P_min
                    self.drawColor2D(scalar_val, frame_title, cur_max, cur_min)
                elif start_end_label=='jacobi_loop':
                    frame_title += ' jacobi loop: ' + str(frame_i)
                    self.drawColor2D(scalar_val, frame_title, P_max, P_min)
                elif start_end_label=='jacobi_senddown':
                    frame_title += ' sent down'
                    self.drawColor2D(scalar_val, frame_title, P_max, P_min)
                elif start_end_label=='jacobi_frame':
                    frame_title += ' after jacobi for loop: '+str(numLoops)
                    self.drawColor2D(scalar_val,frame_title, P_max, P_min)
            else:
                frame_title = label + ' ' + start_end_label
                DivW_dat_i=0
                if label == 'DivW':
                    scalar_val = data_stack[i]
                    DivW_dat_i=i
                    self.drawColor2D(scalar_val, frame_title, DivW_max, DivW_min)
                elif label == 'gradP':
                    scalar_val = data_stack[i]
                    self.drawColor2D(scalar_val, frame_title, gradP_max, gradP_min)
                elif label == 'LapP':
                    scalar_val = data_stack[i]
                    self.drawColor2D(scalar_val, frame_title, DivW_max, DivW_min)
                    DivW_dat = data_stack[DivW_dat_i]
                    if(len(scalar_val)==len(DivW_dat)):
                        frame_title += ' diff with DivW'
                        error_jac = [val1 - val2 for val1, val2 in zip(DivW_dat, scalar_val)]
                        self.drawColor2D(error_jac, frame_title, DivW_max, DivW_min)

        return 0
    
    def drawTestCheck(self, header_stack, data_stack):
        W_min, W_max, DivW_min, DivW_max, P_min, P_max, gradP_min, gradP_max = self.get_W_DivW_P_GradP_limits(header_stack, data_stack) 
        i_stack_base=-1
        i_stack_height=0
        i_P_final_frame=-1
        divW_width=0
        divW_height=0
        P_width=0
        P_height=0
        for i in range(len(header_stack)):
            header = header_stack[i]
            grid_width = ct.getGridWidth(header)
            grid_height = ct.getGridHeight(header)
            label = ct.header2[ct.getDataLabel(header)]
            start_end_label = ct.header4[ct.getDataStartEndCode(header)]
            jacobi_frame_label = 'jacobi'
            if start_end_label=='jacobi_stack_fill':
                stack_i = ct.getJacobiStackI(header)
                if label == 'DivW':
                    i_stack_base=i
                    i_stack_height=stack_i
                    divW_width=grid_width
                    divW_height = grid_height
            elif start_end_label=='jacobi_frame':
                if label == 'P':
                    i_P_final_frame=i
                    P_width = grid_width
                    P_height = grid_height
        i_stack_height+=1
                    
        testStack = tc.genReducedStack(data_stack[i_stack_base], divW_width, divW_height,i_stack_height)
        delta_x = 1e-3
        LapP = tc.Lap(data_stack[i_P_final_frame], P_width, P_height, delta_x)
        max_stack_index = i_stack_height-1    

        for i in range(len(header_stack)):
            header = header_stack[i]
            grid_width = ct.getGridWidth(header)
            grid_height = ct.getGridHeight(header)
            self.setGridWH(grid_width, grid_height)
            label = ct.header2[ct.getDataLabel(header)]
            start_end_label = ct.header4[ct.getDataStartEndCode(header)]
            jacobi_frame_label = 'jacobi'
            if jacobi_frame_label in start_end_label:
                stack_i = ct.getJacobiStackI(header)
                frame_title = label + "  stack: " +str(stack_i)
                scalar_val = data_stack[i]
                if start_end_label=='jacobi_stack_fill':
                    frame_title += ' fill'
                    cur_max=0.0
                    cur_min=0.0
                    test_stack_i = max_stack_index-stack_i
                    test_stack_data = testStack[test_stack_i]
                    if label == 'DivW':
                        cur_max = DivW_max
                        cur_min = DivW_min
                        self.drawColor2D(scalar_val, frame_title, cur_max, cur_min)
                        test_frame_title = frame_title+' test'
                        self.drawColor2D(test_stack_data, test_frame_title, cur_max, cur_min)
            else:
                if label == 'LapP':
                    frame_title = label
                    scalar_val = data_stack[i]
                    self.drawColor2D(scalar_val, frame_title, DivW_max, DivW_min)
                    test_frame_title = "test LapP"
                    self.drawColor2D(LapP,test_frame_title,DivW_max, DivW_min)
                    if(len(scalar_val)==len(LapP)):
                        test_frame_title += ' diff LapP'
                        error_jac = [val1 - val2 for val1, val2 in zip(LapP, scalar_val)]
                        self.drawColor2D(error_jac, test_frame_title, DivW_max, DivW_min)
        
            

class Test:

    def __init__(self, filename = '../Dat/frames.dat'):#'../Dat/test.dat'):#'../Dat/frames.dat'):
        self.filename = filename
        self.C_Trans = ct.CTrans(filename)
        self.Draw = drawFrame(0)

    def testRun(self):
        dat, file_offset = self.C_Trans.readDdat(4,0)
        print(dat)

    def Run(self):
        len_returned = 2
        frame_cnt=0
        while len_returned>1:
            header_stack, data_stack, new_file_offset = self.C_Trans.readFrame()
            len_returned = len(header_stack)
            self.Draw.updateFrameIndex(frame_cnt)
            frame_cnt += 1
            if(len_returned>=1):
                testtt=1
                self.Draw.draw(header_stack, data_stack)
                #self.Draw.drawJacobi(header_stack, data_stack)
                #self.Draw.drawTestCheck(header_stack, data_stack)
            
testInst = Test()
testInst.Run()
print("finished")