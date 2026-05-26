#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#ifndef BASE_H
#include "Base.h"
#endif


__device__ int getIndex(const int i, const int j, const int grid_width) {
    return j * grid_width + i;
}
__device__ int getGoodIndex(const int i, const int j, const int grid_width, const int grid_height) {
    if (i < 0 || j < 0)
        return -1;
    if (i >= grid_width || j >= grid_height)
        return -1;
    return getIndex(i, j, grid_width);
}

__global__ void copy_memory(double* copy, const double* orig) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    copy[i] = orig[i];
}
__global__ void add_values(double* in_out, const double* in) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
	in_out[i] += in[i];
}

/*used for test*/
__global__ void derivative_Core(
    double* dx_out,
    double* dy_out,
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
    dx_out[center_index] = dW_dx;
    dy_out[center_index] = dW_dy;
}
/*used for test*/
__global__ void gradient_core(
    double* Vx,
    double* Vy,
    const double* p,
    double inv_2delta_x,
    const int grid_width,
    const int grid_height)
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
        px = 0.0;
    if (j < (grid_height - 1) && j>0) {
        int index_jplus = getIndex(i, j + 1, grid_width);
        int index_jminus = getIndex(i, j - 1, grid_width);
        py = (p[index_jplus] - p[index_jminus]) * inv_2delta_x;
    }
    else
        py = 0.0;
    Vx[field_index] = px;
    Vy[field_index] = py;
}