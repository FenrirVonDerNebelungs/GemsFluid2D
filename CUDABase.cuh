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
__device__ bool inInteriorSlab_ij(const int i, const int j, const int grid_width, const int grid_height) {
    if (i <= 0 || i >= (grid_width - 1))
        return false;
    if (j <= 0 || j >= (grid_height - 1))
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
    const double Q_22)
{
    double est_val = dist_to_x2 * dist_to_y2 * Q_11 + dist_to_x1 * dist_to_y2 * Q_21 + dist_to_x2 * dist_to_y1 * Q_12 + dist_to_x1 * dist_to_y1 * Q_22;
    return est_val;
}
__device__ bool bilinearAprox_Core(
    double2* pU_out,
    const double* Ux_in,
    const double* Uy_in,
    const double2 relPos,
    const int i_pos,
	const int j_pos,
    const int grid_width, 
	const int grid_height
) {
	s_four_corners relCorners = getRelFourCorners(relPos);
	s_four_corners corners = addCornersTo_ij(i_pos, j_pos, relCorners);

    double dist_to_x1 = relPos.x - (double)relCorners.i1;
    double dist_to_x2 = (double)relCorners.i2 - relPos.x;
    double dist_to_y1 = relPos.y - (double)relCorners.j1;
    double dist_to_y2 = (double)relCorners.j2 - relPos.y;

    double Ux_11 = 0.0, Ux_21 = 0.0, Ux_12 = 0.0, Ux_22 = 0.0;
    double Uy_11 = 0.0, Uy_21 = 0.0, Uy_12 = 0.0, Uy_22 = 0.0;
    if (inRange_corners(corners, grid_width, grid_height)) {
        /* frac{1}{(x_2 - x_1)(y_2-y_1)} ((x_2-x)(y_2-y)Q_{11} + (x-x_1)(y_2-y)Q_{21} + (x_2-x)(y-y_1)Q_{12} + (x-x_1)(y-y_1)Q_{22} */
        /* (y_2-y_1)=(x_2-x_1)=delta_x*/

        int index_11 = getIndex(corners.i1, corners.j1, grid_width);
        int index_12 = getIndex(corners.i1, corners.j2, grid_width);
        int index_21 = getIndex(corners.i2, corners.j1, grid_width);
        int index_22 = getIndex(corners.i2, corners.j2, grid_width);

        Ux_11 = Ux_in[index_11];
        Ux_21 = Ux_in[index_21];
        Ux_12 = Ux_in[index_12];
        Ux_22 = Ux_in[index_22];

        Uy_11 = Uy_in[index_11];
        Uy_21 = Uy_in[index_21];
        Uy_12 = Uy_in[index_12];
        Uy_22 = Uy_in[index_22];
    }
    else {
        int index_11 = getGoodIndex(corners.i1, corners.j1, grid_width, grid_height);
        int index_12 = getGoodIndex(corners.i1, corners.j2, grid_width, grid_height);
        int index_21 = getGoodIndex(corners.i2, corners.j1, grid_width, grid_height);
        int index_22 = getGoodIndex(corners.i2, corners.j2, grid_width, grid_height);
        if (index_11 >= 0) {
            Ux_11 = Ux_in[index_11];
            Uy_11 = Uy_in[index_11];
        }
        if (index_12 >= 0) {
            Ux_12 = Ux_in[index_12];
            Uy_12 = Uy_in[index_12];
        }
        if (index_21 >= 0) {
            Ux_21 = Ux_in[index_21];
            Uy_21 = Uy_in[index_21];
        }
        if (index_22 >= 0) {
            Ux_22 = Ux_in[index_22];
            Uy_22 = Uy_in[index_22];
        }
    }
    /*find values in middle of grid where backtrace remains on grid*/
    double estimated_Ux = bilinearAprox_Quad(
        dist_to_x1,
        dist_to_x2,
        dist_to_y1,
        dist_to_y2,
        Ux_11,
        Ux_21,
        Ux_12,
        Ux_22);
    double estimated_Uy = bilinearAprox_Quad(
        dist_to_x1,
        dist_to_x2,
        dist_to_y1,
        dist_to_y2,
        Uy_11,
        Uy_21,
        Uy_12,
        Uy_22);

    pU_out->x = estimated_Ux;
    pU_out->y = estimated_Uy;
    return true;
}
__global__ void bilinearAprox_scaledFrame_Core( /*i and j block indexes run over smaller frame, pre-scaled*/
    double* Ux_out/*blown up*/,
    double* Uy_out/*blown up*/,
    const double* Ux_in,
    const double* Uy_in,
    int grid_width,
    int grid_height,
    int scale_factor) 
{
    if(scale_factor <= 1)
		return;
	double2 U_out = make_double2(0.0, 0.0);
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
	int half_square = scale_factor / 2;
	int out_grid_width = grid_width * scale_factor;
	int out_grid_height = grid_height * scale_factor;
    double2 relPos = make_double2(0.0, 0.0);
    for(int square_j=0; square_j < scale_factor; square_j++) {
		relPos.y = (double)(square_j-half_square) / (double)scale_factor;
        int j_out = j * scale_factor - half_square + square_j;
        if(j_out < 0 || j_out >= out_grid_height)
			continue;
        for(int square_i=0; square_i < scale_factor; square_i++) {
			relPos.x = (double)(square_i - half_square) / (double)scale_factor;
            int i_out = i * scale_factor - half_square + square_i;
            if (i_out < 0 || i_out >= out_grid_width)
                continue;
            int out_index = getIndex(i_out, j_out, out_grid_width);
			bilinearAprox_Core(&U_out, Ux_in, Uy_in, relPos, i, j, grid_width, grid_height);
			Ux_out[out_index] = U_out.x;
			Uy_out[out_index] = U_out.y;
        }
	}
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
