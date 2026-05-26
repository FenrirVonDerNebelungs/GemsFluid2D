#include "Jacobi.cuh"
__device__ bool inRange_corners(const int4* corners, int grid_width, int grid_height) {
    if (corners->x == corners->z || corners->y == corners->w)
        return false;
    return corners->x >= 0 && corners->z < grid_width && corners->y >= 0 && corners->w < grid_height;
}
__device__ void add(const int4* corners1, const int4* corners2, int4* result) {
    result->x = corners1->x + corners2->x;
    result->z = corners1->z + corners2->z;
    result->y = corners1->y + corners2->y;
    result->w = corners1->w + corners2->w;
}
__device__ void addCornersTo_ij(int i, int j, const int4* relCorners, int4* result) {
    int4 indexes_to_add = make_int4(i, j, i, j); /* x, y, z, w*/
    add(&indexes_to_add, relCorners, result);
}
__device__ void getRelFourCorners(
    const double2& relPos,
    int4* result
) {
    result->x = (int)std::floor(relPos.x);
    result->z = (int)std::ceil(relPos.x);
    result->y = (int)std::floor(relPos.y);
    result->w = (int)std::ceil(relPos.y);
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

__device__ double bilinearAprox_Core(
    int Q_index,
    const double* Q,
    const int* index_11,
    const int* index_12,
    const int* index_21,
    const int* index_22,
    const double* dist_to_x1,
    const double* dist_to_x2,
    const double* dist_to_y1,
    const double* dist_to_y2
) {
    /* frac{1}{(x_2 - x_1)(y_2-y_1)} ((x_2-x)(y_2-y)Q_{11} + (x-x_1)(y_2-y)Q_{21} + (x_2-x)(y-y_1)Q_{12} + (x-x_1)(y-y_1)Q_{22} */
    /* (y_2-y_1)=(x_2-x_1)=delta_x*/
    double Q_11 = 0.0, Q_21 = 0.0, Q_12 = 0.0, Q_22 = 0.0;
    if (index_11[Q_index] >= 0)
        Q_11 = Q[index_11[Q_index]];
    if (index_21[Q_index] >= 0)
        Q_21 = Q[index_21[Q_index]];
    if (index_12[Q_index] >= 0)
        Q_12 = Q[index_12[Q_index]];
    if (index_22[Q_index] >= 0)
        Q_22 = Q[index_22[Q_index]];
    double estimated_Q = bilinearAprox_Quad(
        dist_to_x1[Q_index],
        dist_to_x2[Q_index],
        dist_to_y1[Q_index],
        dist_to_y2[Q_index],
        Q_11,
        Q_21,
        Q_12,
        Q_22);
    return estimated_Q;
}
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

__global__ void advection_backtrace_Core(
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
__global__ void advection_backtrace_test_Core( /*test function returns backtraced cell center locations based on cell current velocity*/
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