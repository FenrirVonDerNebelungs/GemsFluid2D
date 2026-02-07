#include "CUDABase.cuh"
/* x^{k+1}_{i,j} = frac{x^k_{i-1,j} + x^k_{i+1,j} + x^k_{i,j-1} + x^k_{i,j+1} + \alpha b_{i,j}}{\beta}
* performs jacobi iteration to solve the poisson equation
*/
__global__ void jacobi(
    double* frame_out,
    double* frame_in,
    const double* b_in,
    const double alpha,
    const double rbeta,
    const int grid_width,
    const int grid_height)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (inInteriorSlab_ij(i,j,grid_width, grid_height) ){/*not the most efficent but clear*/
        int center_index = getIndex(i, j, grid_width);
        double xL = frame_in[getIndex(i - 1, j, grid_width)];
        double xR = frame_in[getIndex(i + 1, j, grid_width)];
        double xB = frame_in[getIndex(i, j - 1, grid_width)];
        double xT = frame_in[getIndex(i, j + 1, grid_width)];
        double b = b_in[center_index];
        double xNew = (xL + xR + xB + xT + alpha * b) * rbeta;
        frame_out[center_index] = xNew;
    }
}
/*solve boundary conditions*/
/*solves top bottom, equations of the form
* x_{i, 0} + (-sign)*x_{i,1} = 0
*/
__global__ void boundary_topbottom(double* frame, double sign, const int grid_width, const int grid_height) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    double xB_1 = frame[getIndex(i, 1, grid_width)];
    double xB_0 = xB_1 * sign;
    int out_index = getIndex(i, 0, grid_width);
    frame[out_index] = xB_0;
    double xT_1 = frame[getIndex(i, grid_height - 2, grid_width)];
    double xT_0 = xT_1 * sign;
    out_index = getIndex(i, grid_height - 1, grid_width);
    frame[out_index] = xT_0;
}
/*solves right left equations of the form
* x_{0,j} + (-sign)*x_{1,j} = 0
*/
__global__ void boundary_rightleft(double* frame, double sign, const int grid_width, const int grid_height) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    double xL_1 = frame[getIndex(1, j, grid_width)];
    double xL_0 = xL_1 * sign;
    int out_index = getIndex(0, j, grid_width);
    frame[out_index] = xL_0;
    double xR_1 = frame[getIndex(grid_width - 2, j, grid_width)];
    double xR_0 = xR_1 * sign;
    out_index = getIndex(grid_width - 1, j, grid_width);
    frame[out_index] = xR_0;
}