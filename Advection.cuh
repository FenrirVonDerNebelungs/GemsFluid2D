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
    bilinearAprox_Core(&U_out, Ux_in, Uy_in, relPos, i, j, grid_width, grid_height);
    Ux_out[center_index] = U_out.x;
    Uy_out[center_index] = U_out.y;
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

    s_four_corners relCorners = getRelFourCorners(relPos);
    s_four_corners corners = addCornersTo_ij(i, j, relCorners);
    x_out[center_index] = relPos.x;
    y_out[center_index] = relPos.y;
}