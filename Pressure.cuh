#include "CUDABase.cuh"

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
        WB = Wx_in[index_jminus];
    }
    if (j < j_max) {
        int index_jplus = getIndex(i, j + 1, grid_width);
        WT = Wx_in[index_jplus];
    }
    int center_index = getIndex(i, j, grid_width);
    double dW_dx = inv_2delta_x * (WR - WL);
    double dW_dy = inv_2delta_x * (WT - WB);
    frame_out[center_index] = dW_dx + dW_dy;
}
__device__ void fix_boundary_U(double* Ux, double* Uy, const int field_index, const int i, const int j, const int grid_width, const int grid_height) {
    if (i <= 0 || i >= (grid_width - 1))
        Ux[field_index] = 0.0;
    if (j <= 0 || j >= (grid_height - 1))
        Uy[field_index] = 0.0;
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
    int index_iplus = (i < (grid_width - 1)) ? getIndex(i + 1, j, grid_width) : 0;
    int index_iminus = (i > 0) ? getIndex(i - 1, j, grid_width) : 0;
    int index_jplus = (j < (grid_height - 1)) ? getIndex(i, j + 1, grid_width) : 0;
    int index_jminus = (j > 0) ? getIndex(i, j - 1, grid_width) : 0;
    px = (p[index_iplus] - p[index_iminus]) * inv_2delta_x;
    py = (p[index_jplus] - p[index_jminus]) * inv_2delta_x;
    Ux[field_index] = Wx[field_index] - px;
    Uy[field_index] = Wy[field_index] - py;
    //fix_boundary_U(Ux, Uy, field_index, i, j, grid_width, grid_height);/* Ux and Uy should already be 0 along the corresponding boundaries, but this corrects for rounding errors and the corners*/
}