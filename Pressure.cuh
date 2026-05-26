#include "Forces.cuh"
#define PRESSURE_CUH_EDGE_RED_LEN 16

__device__ double edge_pressure_average;

__global__ void divergence_Core(
    double* frame_out,
    const double* Wx_in,
    const double* Wy_in,
    double inv_2delta_x,
    const int grid_width,
    const int grid_height) /*this is used to compute div of velocity field, so assume dirichlet b.c. w_{-1,j} = 0 */
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    double WL = 0.0, WR = 0.0;
    double WT = 0.0, WB = 0.0;
    int i_max = grid_width - 1;
    int j_max = grid_height - 1;
    if (i > 0) {
        int index_iminus = getIndex(i - 1, j, grid_width);
        WL = Wx_in[index_iminus];
    }
    if (i < i_max) {
        int index_iplus = getIndex(i + 1, j, grid_width);
        WR = Wx_in[index_iplus];
    }
    if (j > 0) {
        int index_jminus = getIndex(i, j - 1, grid_width);
        WB = Wy_in[index_jminus];
    }
    if (j < j_max) {
        int index_jplus = getIndex(i, j + 1, grid_width);
        WT = Wy_in[index_jplus];
    }
    int center_index = getIndex(i, j, grid_width);
    double dW_dx = inv_2delta_x * (WR - WL);
    double dW_dy = inv_2delta_x * (WT - WB);
    frame_out[center_index] = dW_dx + dW_dy;
}
__global__ void subtractGradient(
    double* Ux,
    double* Uy,
    const double* Wx,
    const double* Wy,
    const double* p,
    double inv_2delta_x,
    const int grid_width,
    const int grid_height) /* fixes w->u everwhere except at the four corners */
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int field_index = getIndex(i, j, grid_width);
    double px = 0.0;
    double py = 0.0;
    if (i < (grid_width - 1) && i>0) {
        int index_iplus = getIndex(i + 1, j, grid_width);
        int index_iminus = getIndex(i - 1, j, grid_width);
        px = (p[index_iplus] - p[index_iminus]) * inv_2delta_x;
    }
    else
        px = Wx[field_index];
    if (j < (grid_height - 1) && j>0) {
        int index_jplus = getIndex(i, j + 1, grid_width);
        int index_jminus = getIndex(i, j - 1, grid_width);
        py = (p[index_jplus] - p[index_jminus]) * inv_2delta_x;
    }
    else
        py = Wy[field_index];

    Ux[field_index] = Wx[field_index] - px;
    Uy[field_index] = Wy[field_index] - py;
}
__global__ void extract_edges(double* X_1D, const double* X_2D, int width, int height) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j_max = (height - 1) * width;
    int i_max = (width - 1);
    int index2D_j0 = i;
    int index2D_jmax = j_max + i;
    int index2D_i0 = i * width;/* j*width */
    int index2D_imax = index2D_i0 + i_max; /* j*width + i_max */
    
    int index1D_j0 = i;
    int index1D_jmax = i + width;
    int index1D_i0 = 2 * width + i;
    int index1D_imax = 2 * width + height + i;
    
    X_1D[index1D_j0] = X_2D[index2D_j0];
    X_1D[index1D_jmax] = X_2D[index2D_jmax];
    X_1D[index1D_i0] = X_2D[index2D_i0];
    X_1D[index1D_imax] = X_2D[index2D_imax];
    if (i == 0)
        X_1D[index1D_i0] = 0.0;
    else if (i == i_max)
        X_1D[index1D_imax] = 0.0;
}
__global__ void average_array(double* X_out, const double* X_in, int reduction_factor) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int start_in_index = i * reduction_factor;
    double sum_val = 0.0;
    double red_factor = static_cast<double>(reduction_factor);
    for (int red_i = 0; red_i < reduction_factor; red_i++) {
        sum_val += X_in[start_in_index + red_i];
    }
    X_out[i] = sum_val / red_factor;
}
__global__ void final_average_array(const double* X_in, int X_in_len) {
    //int i = blockIdx.x * blockDim.x + threadIdx.x;/*only index i=0*/
    double sum_val = 0.0;
    double num_in_ave = 0.0;
    for (int red_i = 0; red_i < X_in_len; red_i++) {
        sum_val += X_in[red_i];
        num_in_ave += 1.0;
    }
   edge_pressure_average = sum_val/num_in_ave;
}
__global__ void reset_net_pressure_Core(double* X_in_out) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    X_in_out[i] -= edge_pressure_average;
}