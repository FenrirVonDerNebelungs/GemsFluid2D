#include "Pressure.cuh"
/* x^{k+1}_{i,j} = frac{x^k_{i-1,j} + x^k_{i+1,j} + x^k_{i,j-1} + x^k_{i,j+1} + \alpha b_{i,j}}{\beta}
* performs jacobi iteration to solve the poisson equation
*/

__global__ void jacobi(
    double* frame_out,
    const double* frame_in,
    const double* b_in,
    const double alpha,
    const double rbeta,
    const int grid_width,
    const int grid_height) /* without correction fullfills dirichlet boundary conditions with x=0 for values off grid*/
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int center_index = getIndex(i, j, grid_width);
    double xL = 0.0, xR = 0.0;
    double xT = 0.0, xB = 0.0;

    if(i>0)
        xL = frame_in[getIndex(i - 1, j, grid_width)];
    if(i<(grid_width-1))
        xR = frame_in[getIndex(i + 1, j, grid_width)];
    if(j>0)
        xB = frame_in[getIndex(i, j - 1, grid_width)];
    if(j<(grid_height-1))
        xT = frame_in[getIndex(i, j + 1, grid_width)];
    
    double b = b_in[center_index];
    double xNew = (xL + xR + xB + xT + alpha * b) * rbeta;
    frame_out[center_index] = xNew;
}
__device__ void jacobi_boundary_pressure_setVal(
    double* frame_out,
    const double* b_in,
    int i,
    int j,
    double xL,
    double xR,
    double xB,
    double xT,
    double g,
    double alpha,
    double rbeta,
    int grid_width)
{
    int current_index = getIndex(i, j, grid_width);
    double b = b_in[current_index];
    double xNew = (xL + xR + xB + xT - g + alpha * b) * rbeta;
    frame_out[current_index] = xNew;
}
__device__ void jacobi_boundary_pressure_doSide(
    double* frame_out,
    const double* frame_in,
    const double* W_in,
    const double * b_in,
    int i,
    int iL,
    int iR,
    int j,
    int jT,
    int jB,
    double alpha,
    double rbeta,
    double delta_x,
    int grid_width,
    int grid_height)
{
    //double alpha = -delta_x * delta_x;
    int i_max = grid_width - 1;
    int j_max = grid_height - 1;
    double xL = 0.0, xR = 0.0;
    double xT = 0.0, xB = 0.0;
    double g = 0.0;
    int g_index = 0;
    double g_mult_constant = 2.0 * delta_x;
    if(jB>=0)
        xB = frame_in[getIndex(i, jB, grid_width)];
    if(jT>=0)
        xT = frame_in[getIndex(i, jT, grid_width)];
    if(iL>=0)
        xL = frame_in[getIndex(iL, j, grid_width)];
    if(iR>=0)
        xR = frame_in[getIndex(iR, j, grid_width)];

    if(jB<0){
        g_index = getIndex(i, 0, grid_width);
        xB = xT;
    }
    if(jT<0){
        g_index = getIndex(i, j_max, grid_width);
        g_mult_constant = -g_mult_constant;
        xT = xB;
    }
    if(iL<0){
        g_index = getIndex(0, j, grid_width);
        xL = xR;
    }
    if(iR<0){
        g_index = getIndex(i_max, j, grid_width);
        g_mult_constant = -g_mult_constant;
        xR = xL;
    }
    g = g_mult_constant * W_in[g_index];
    jacobi_boundary_pressure_setVal(frame_out, b_in, i, j, xL, xR, xB, xT, g, alpha, rbeta, grid_width);
}
__global__ void jacobi_boundary_pressure(
    double* frame_out,
    const double* Wx_in,
    const double* Wy_in,
    const double* frame_in /* p */,
    const double* b_in/* div * w */,
    const double alpha,
    const double rbeta,
    const double delta_x,
    const int grid_width,
    const int grid_height)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    int i = 0, j = 0;
    int iL = 0, iR = 0;
    int jB = 0, jT = 0;
    /** use k as j and run along x=0 and x=max **/
    if (k > 0 && k < (grid_height - 1)) {
        j = k;
        jB = j - 1;
        jT = j + 1;
        /* go along i=1 */
        i = 0;
        iR = i + 1;
        iL = -1;
        jacobi_boundary_pressure_doSide(frame_out, frame_in, Wx_in, b_in, i, iL, iR, j, jT, jB, alpha, rbeta, delta_x, grid_width, grid_height);
        /* go along i=i_max-1 */
        i = grid_width - 1;
        iR = -1;
        iL = i - 1;
        jacobi_boundary_pressure_doSide(frame_out, frame_in, Wx_in, b_in, i, iL, iR, j, jT, jB, alpha, rbeta, delta_x, grid_width, grid_height);
    //}
    /** use k as i and run along y=0 and y=max **/
    //if (k > 0 && k < (grid_width - 1)) {
        i = k;
        iL = i - 1;
        iR = i + 1;
        /* go along j=1 */
        j = 0;
        jT = j + 1;
        jB = -1;
        jacobi_boundary_pressure_doSide(frame_out, frame_in, Wy_in, b_in, i, iL, iR, j, jT, jB, alpha, rbeta, delta_x, grid_width, grid_height);
        /* go along j=j_max-1 */
        j = grid_height - 1;
        jT = -1;
        jB = j - 1;
        jacobi_boundary_pressure_doSide(frame_out, frame_in, Wy_in, b_in, i, iL, iR, j, jT, jB, alpha, rbeta, delta_x, grid_width, grid_height);
    }
}
