#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#ifndef BASE_H
#include "Base.h"
#endif


__device__ int getIndex(const int i, const int j, const int grid_width) {
    return j * grid_width + i;
}
__device__ bool inInteriorSlab_ij(const int i, const int j, const int grid_width, const int grid_height) {
    if (i < 1 || i >= (grid_width - 1))
        return false;
    if (j < 1 || j >= (grid_height - 1))
        return false;
    return true;
}
__device__ bool inRange_corners(const s_four_corners& corners, int grid_width, int grid_height) {
    return corners.i1 >= 0 && corners.i2 < grid_width && corners.j1 >= 0 && corners.j2 < grid_height;
}
__device__ s_four_corners add(const s_four_corners& corners1, const s_four_corners& corners2) {
    s_four_corners corners = { corners1.i1 + corners2.i1, corners1.i2 + corners2.i2, corners1.j1 + corners2.j1, corners1.j2 + corners2.j2 };
    return corners;
}
__device__ s_four_corners addCornersTo_ij(int i, int j, const s_four_corners& relCorners) {
    s_four_corners indexes_to_add = { i, i, j, j };
    s_four_corners corners = add(indexes_to_add, relCorners);
    return corners;
}
__device__ s_four_corners getRelFourCorners(
    const int i,
    const int j,
    const double2& relPos
) {
    s_four_corners corners;
    corners.i1 = (int)std::floor(relPos.x);
    corners.i2 = (int)std::ceil(relPos.x);
    corners.j1 = (int)std::floor(relPos.y);
    corners.j2 = (int)std::ceil(relPos.y);
    return corners;
}
__device__ double bilinearAprox_Quad(
    const double dist_to_x1, 
    const double dist_to_x2, 
    const double dist_to_y1, 
    const double dist_to_y2,
    const double Q_11,
    const double Q_21,
    const double Q_12,
    const double Q_22,
    const double delta_x_sqr_inv)
{
    double est_val = dist_to_x2 * dist_to_y2 * Q_11 + dist_to_x1 * dist_to_y2 * Q_21 + dist_to_x2 * dist_to_y1 * Q_12 + dist_to_x1 * dist_to_y1 * Q_22;
    est_val *= delta_x_sqr_inv;
    return est_val;
}
__device__ double2 bilinearAprox_Core(
    const double* Ux_in,
    const double* Uy_in,
    const double2 relPos,
    const s_four_corners& relCorners,
    const s_four_corners& corners,
    const double delta_x,
    const int grid_width
) {
    /* frac{1}{(x_2 - x_1)(y_2-y_1)} ((x_2-x)(y_2-y)Q_{11} + (x-x_1)(y_2-y)Q_{21} + (x_2-x)(y-y_1)Q_{12} + (x-x_1)(y-y_1)Q_{22} */
    /* (y_2-y_1)=(x_2-x_1)=delta_x*/
    double dist_to_x1 = relPos.x - (double)relCorners.i1;
    double dist_to_x2 = (double)relCorners.i2 - relPos.x;
    double dist_to_y1 = relPos.y - (double)relCorners.j1;
    double dist_to_y2 = (double)relCorners.j2 - relPos.y;
    s_four_corners goodcorners = corners;
    int index_11 = getIndex(goodcorners.i1, goodcorners.j1, grid_width);
    int index_12 = getIndex(goodcorners.i1, goodcorners.j2, grid_width);
    int index_21 = getIndex(goodcorners.i2, goodcorners.j1, grid_width);
    int index_22 = getIndex(goodcorners.i2, goodcorners.j2, grid_width);
    double estimated_Ux, estimated_Uy;

    /*find values in middle of grid where backtrace remains on grid*/
    double delta_x_sqr_inv = 1.0 / (delta_x * delta_x);
    estimated_Ux = bilinearAprox_Quad(
        dist_to_x1,
        dist_to_x2,
        dist_to_y1,
        dist_to_y2,
        Ux_in[index_11],
        Ux_in[index_21],
        Ux_in[index_12],
        Ux_in[index_22],
        delta_x_sqr_inv);
    estimated_Uy = bilinearAprox_Quad(
        dist_to_x1,
        dist_to_x2,
        dist_to_y1,
        dist_to_y2,
        Uy_in[index_11],
        Uy_in[index_21],
        Ux_in[index_12],
        Ux_in[index_22],
        delta_x_sqr_inv);

    double2 U_out = make_double2(estimated_Ux, estimated_Uy);
    return U_out;
}
__global__ void copy_memory(double* copy, const double* orig) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    copy[i] = orig[i];
}
__global__ void applyForce_Core(
    double* Ux_out,
    double* Uy_out,
    const double* Ux_in,
    const double* Uy_in,
    const double2 center/*in i,j*/,
    const double2 F_c,
    double delta_t,
    double inv_rsqrd,
    const int grid_width) 
{/*currently runs over entire frame, this could be reduced*/
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int center_index = getIndex(i, j, grid_width);
    /* frac{(x-x_p)^2 + (y-y_p)^2}{r}*/
    double x_dist = (double)i - center.x;
    double y_dist = (double)j - center.y;
    double arg = -(x_dist * x_dist + y_dist * y_dist) * inv_rsqrd;
    if (arg > -4.6) {
        double force_x = F_c.x * exp(arg);
        double force_y = F_c.y * exp(arg);
        double delta_ux = delta_t * force_x;
        double delta_uy = delta_t * force_y;
        Ux_out[center_index] = Ux_in[center_index] + delta_ux;
        Uy_out[center_index] = Uy_in[center_index] + delta_uy;
    }
    else {
        Ux_out[center_index] = Ux_in[center_index];
        Uy_out[center_index] = Uy_in[center_index];
    }
}
__global__ void divergence_Core(
    double* frame_out, 
    const double* Wx_in, 
    const double* Wy_in, 
    double inv_2delta_x, 
    const int grid_width, 
    const int grid_height) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= 1 && i < (grid_width - 1) && j >= 1 && j < (grid_height - 1)) {/*not the most efficent but clear*/
        int center_index = getIndex(i, j, grid_width);
        int index_iplus = getIndex(i + 1, j, grid_width);
        int index_iminus = getIndex(i - 1, j, grid_width);
        int index_jplus = getIndex(i, j + 1, grid_width);
        int index_jminus = getIndex(i, j - 1, grid_width);
        double dW_dx = inv_2delta_x * (Wx_in[index_iplus] - Wx_in[index_iminus]);
        double dW_dy = inv_2delta_x * (Wy_in[index_jplus] - Wy_in[index_jminus]);
        frame_out[center_index] = dW_dx + dW_dy;
    }
}
__global__ void subtractGradient(
    double* Ux, 
    double* Uy, 
    const double* Wx, 
    const double* Wy, 
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
    if (i >= 1 && i < (grid_width - 1)) {
        int index_iplus = getIndex(i + 1, j, grid_width);
        int index_iminus = getIndex(i - 1, j, grid_width);
        px = (p[index_iplus] - p[index_iminus]) * inv_2delta_x;
    }
    if (j >= 1 && j < (grid_height - 1)) {
        int index_jplus = getIndex(i, j + 1, grid_width);
        int index_jminus = getIndex(i, j - 1, grid_width);
        py = (p[index_jplus] - p[index_jminus]) * inv_2delta_x;
    }
    Ux[field_index] = Wx[field_index] - px;
    Uy[field_index] = Wy[field_index] - py;
}