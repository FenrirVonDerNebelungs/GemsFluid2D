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