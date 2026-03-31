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
    def __init__(self):
        self.grid_wh = g_grid_width
        self.grid_exp_wh = g_grid_width*g_blow_factor

    def convertTo2D(self, raw_array):
        V_array = np.zeros((g_grid_width, g_grid_width))
        for j in range(g_grid_width):
            for i in range(g_grid_width):
                V_array[i,j]=raw_array[j*g_grid_width + i]
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

    def draw_scalar(self, p, graph_title):
        x,y = np.meshgrid(np.arrange(g_delta_x, g_grid_width*g_delta_x, g_delta_x), np.arrange(g_delta_x, g_grid_width*g_delta_x, g_delta_x))
        p_2D = self.convertTo2D(p)
        fig=plt.figure(figsize=(7,7))
        ax=plt.axes(projection='3d')
        surf = ax.plot_surface(x,y,p_2D,cmap='coolwarm')
        ax.set_title(graph_title)
        plt.show()

    def draw(header_stack, data_stack):
        return 0

class Test:

    def __init__(self, filename = '../Dat/frames.dat'):
        grid_stream_len = g_grid_width*g_grid_width
        grid_exp_stream_len = grid_stream_len*g_blow_factor*g_blow_factor
        total_grid_stream_len = grid_stream_len*(11+g_num_jacobi_frames)+grid_exp_stream_len*2
        total_num_stream_headers = (11+g_num_jacobi_frames)+2
        header_len_bytes = ct.header_len
        total_grid_stream_len_bytes = total_grid_stream_len*8
        total_headers_len_bytes = header_len_bytes*total_num_stream_headers
        self.framesize = total_grid_stream_len_bytes + total_headers_len_bytes
        self.filename = filename
        