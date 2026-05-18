#include "Jacobi.cuh"
__device__ double2 getRelBackTracedPosition(
    double2 U /*result is in i,j coordinates where delta_x between cells is 1 */,
    const double delta_t,
    const double delta_x) 
{
    double Ux_back = -U.x * delta_t / delta_x;
    double Uy_back = -U.y * delta_t / delta_x;
    double2 relPos = make_double2(Ux_back, Uy_back); 
    return relPos;
}

__global__ void advection_Core(
    double* Ux_out,
    double* Uy_out,
    double* dye_out,
    const double* Ux_in,
    const double* Uy_in,
    const double* dye_in,
    const int* index_11,
    const int* index_12,
    const int* index_21,
    const int* index_22,
    const double* dist_to_x1,
    const double* dist_to_x2,
    const double* dist_to_y1,
    const double* dist_to_y2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    /* find the interpolated values for U_x, U_y*/
	double est_Ux = bilinearAprox_Core(i, Ux_in, index_11, index_12, index_21, index_22, dist_to_x1, dist_to_x2, dist_to_y1, dist_to_y2);
	double est_Uy = bilinearAprox_Core(i, Uy_in, index_11, index_12, index_21, index_22, dist_to_x1, dist_to_x2, dist_to_y1, dist_to_y2);
	double est_dye = bilinearAprox_Core(i, dye_in, index_11, index_12, index_21, index_22, dist_to_x1, dist_to_x2, dist_to_y1, dist_to_y2);
    Ux_out[i] = est_Ux;
    Uy_out[i] = est_Uy;
    dye_out[i] = est_dye;
}
__global__ void advection_Core_big(
    double* Ux_out, 
    double* Uy_out, 
    const double* Ux_in, 
    const double* Uy_in, 
    double delta_t, 
    double delta_x, 
    int grid_width, 
    int grid_height) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int center_index = getIndex(i, j, grid_width);
    double2 U_in = make_double2(Ux_in[center_index], Uy_in[center_index]);
    double2 relPos = getRelBackTracedPosition(U_in, delta_t, delta_x);
    double2 U_out=make_double2(0.0, 0.0);
    /* find the interpolated values for U_x, U_y*/
    bilinearAprox_Core_big(&U_out, Ux_in, Uy_in, relPos, i, j, grid_width, grid_height);
    Ux_out[center_index] = U_out.x;
    Uy_out[center_index] = U_out.y;
}

__global__ void advection_backtrace_to_indexes(
    int* index_11,
    int* index_12,
    int* index_21,
    int* index_22,
    double* dist_to_x1,
    double* dist_to_x2,
    double* dist_to_y1,
    double* dist_to_y2,
    const double* Ux_in,
    const double* Uy_in,
    double delta_t,
    double delta_x,
    int grid_width,
    int grid_height)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int center_index = getIndex(i, j, grid_width);
    double2 U_in = make_double2(Ux_in[center_index], Uy_in[center_index]);
    double2 relPos = getRelBackTracedPosition(U_in, delta_t, delta_x);
    double2 U_out = make_double2(0.0, 0.0);
    int4 relCorners = make_int4(0, 0, 0, 0);
    getRelFourCorners(relPos, &relCorners);
    int4 corners = make_int4(0, 0, 0, 0);
    addCornersTo_ij(i, j, &relCorners, &corners);

    /*x is i1, z is i2, y is j1, w is j2*/
    dist_to_x1[center_index] = relPos.x - (double)relCorners.x;
    dist_to_x2[center_index] = (double)relCorners.z - relPos.x;
    dist_to_y1[center_index] = relPos.y - (double)relCorners.y;
    dist_to_y2[center_index] = (double)relCorners.w - relPos.y;
    if (inRange_corners(&corners, grid_width, grid_height)) {
        index_11[center_index] = getIndex(corners.x, corners.y, grid_width); /* i1, j1*/
        index_12[center_index] = getIndex(corners.x, corners.w, grid_width); /* i1, j2*/
        index_21[center_index] = getIndex(corners.z, corners.y, grid_width); /* i2, j1*/
        index_22[center_index] = getIndex(corners.z, corners.w, grid_width); /* i2, j2*/
    }
    else {
        index_11[center_index] = getGoodIndex(corners.x, corners.y, grid_width, grid_height);
        index_12[center_index] = getGoodIndex(corners.x, corners.w, grid_width, grid_height);
        index_21[center_index] = getGoodIndex(corners.z, corners.y, grid_width, grid_height);
		index_22[center_index] = getGoodIndex(corners.z, corners.w, grid_width, grid_height);
        if (corners.x == corners.z) {
            dist_to_x1[center_index] = 0.5;
            dist_to_x2[center_index] = 0.5;
        }
        if (corners.y == corners.w) {
            dist_to_y1[center_index] = 0.5;
            dist_to_y2[center_index] = 0.5;
        }
    }
}
__global__ void advection_backtrace_Core( /*test function returns backtraced cell center locations based on cell current velocity*/
    double* x_out,
    double* y_out,
    const double* Ux_in,
    const double* Uy_in,
    double delta_t,
    double delta_x,
    int grid_width,
    int grid_height)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int center_index = getIndex(i, j, grid_width);
    double2 U_in = make_double2(Ux_in[center_index], Uy_in[center_index]);
    double2 relPos = getRelBackTracedPosition(U_in, delta_t, delta_x);

	int4 relCorners = make_int4(0, 0, 0, 0);
    getRelFourCorners(relPos, &relCorners);
	int4 corners = make_int4(0, 0, 0, 0);
    addCornersTo_ij(i, j, &relCorners, &corners);
    x_out[center_index] = relPos.x;
    y_out[center_index] = relPos.y;
}