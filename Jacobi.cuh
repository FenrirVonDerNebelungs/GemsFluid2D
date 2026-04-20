#include "Pressure.cuh"

__device__ int getReducedIndex(const int i, const int j, const int r_grid_width) {
	int reduction_i = floor(i / 2);
	int reduction_j = floor(j / 2);
    return reduction_j * r_grid_width + reduction_i;
}
__device__ int getExpandedIndex(const int r_i, const int r_j, const int offset_i, const int offset_j, const int grid_width) {
    int expansion_i = r_i * 2 +offset_i;
    int expansion_j = r_j * 2+offset_j;
    return expansion_j * grid_width + expansion_i;
}
/*r_X_corner has length of 4*/
__device__ void getGridCorners(double r_X_corner[], const double* r_grid, const int r_i, const int r_j, const int r_grid_width){
	r_X_corner[0] = r_grid[getIndex(r_i, r_j, r_grid_width)];
	r_X_corner[1] = r_grid[getIndex(r_i + 1, r_j, r_grid_width)];
	r_X_corner[2] = r_grid[getIndex(r_i, r_j + 1, r_grid_width)];
	r_X_corner[3] = r_grid[getIndex(r_i + 1, r_j + 1, r_grid_width)];
}
/*r_grid is 1/2 width(height) of grid
* i=2I+1, j=2J+1, for grid cell at center of each reduced grid cell
*/
__device__ void Xgrid_reduction(
    double* r_grid,
    const int r_index,
    const double* grid,
    const int i_center,
    const int i_L,
    const int i_R,
    const int i_B,
    const int i_T) {
    const double center_w = 1.0 / 2.0;
    const double side_w = 1.0 / 4.0;
    const double ave_w = 1.0 / (center_w + 4.0 * side_w);
    double p_ij = grid[i_center];
    double p_side_sum = grid[i_L] + grid[i_R] + grid[i_B] + grid[i_T];
    p_ij *= center_w;
    p_side_sum *= side_w;
    r_grid[r_index] = ave_w * (p_ij + p_side_sum);
}
/*r_grid must be 1/2 width and 1/2 height of grid
* blocks/threads run over REDUCED grid */
__global__ void Xgrid_reduction_x3(
    double* r_grid1,
	double* r_grid2,
	double* r_grid3,
    const int r_grid_width, 
    const int r_grid_height, 
    const double* grid1,
    const double* grid2,
	const double* grid3,
    const int grid_width, 
    const int grid_height)
{
    int r_i = blockIdx.x * blockDim.x + threadIdx.x;
    int r_j = blockIdx.y * blockDim.y + threadIdx.y;
    int r_index = getIndex(r_i, r_j, r_grid_width);
    int i = 2 * r_i + 1;
    int j = 2 * r_j + 1;
    int i_center = getIndex(i, j, grid_width);
    int i_L = getIndex(i - 1, j, grid_width);
    int i_R = getIndex(i + 1, j, grid_width);
    int i_B = getIndex(i, j - 1, grid_width);
    int i_T = getIndex(i, j + 1, grid_width);
	Xgrid_reduction(r_grid1, r_index, grid1, i_center, i_L, i_R, i_B, i_T);
	Xgrid_reduction(r_grid2, r_index, grid2, i_center, i_L, i_R, i_B, i_T);
	Xgrid_reduction(r_grid3, r_index, grid3, i_center, i_L, i_R, i_B, i_T);
}
/*grid must be 2*width and 2*height of r_grid
* blocks/threads run over REDUCED grid */
__global__ void Xgrid_expansion(double* grid, const int grid_width, const int grid_height, const double* r_grid, const int r_grid_width, const int r_grid_height)
{
    int r_i = blockIdx.x * blockDim.x + threadIdx.x;
    int r_j = blockIdx.y * blockDim.y + threadIdx.y;

    /**fill center**/
    int r_index = getIndex(r_i, r_j, r_grid_width);
    int i_center = 2 * r_i + 1;
    int j_center = 2 * r_j + 1;
    grid[getIndex(i_center, j_center, grid_width)] = r_grid[r_index];
    /*fill rest of central blocks*/
    if(inInteriorSlab_ij(r_i, r_j, r_grid_width, r_grid_height)){
        /** fill center of side blocks **/
        int index_R = getIndex(i_center + 1, j_center, grid_width);
        grid[index_R] = 0.5 * (r_grid[r_index] + r_grid[getIndex(r_i + 1, r_j, r_grid_width)]);
        int index_T = getIndex(i_center, j_center + 1, grid_width);
        grid[index_T] = 0.5 * (r_grid[r_index] + r_grid[getIndex(r_i, r_j + 1, r_grid_width)]);
        /** fill corner blocks  filling top right block**/
        int index_RT = getIndex(i_center + 1, j_center + 1, grid_width);
        grid[index_RT] = 0.25 * (r_grid[r_index] + r_grid[getIndex(r_i + 1, r_j, r_grid_width)] + r_grid[getIndex(r_i, r_j + 1, r_grid_width)] + r_grid[getIndex(r_i + 1, r_j + 1, r_grid_width)]);
    }
    else {
        /**fill along left side blocks**/
        if (r_i == 0) {
            /*** fill right and left center side blocks ***/
            int index_R = getIndex(i_center + 1, j_center, grid_width);
            grid[index_R] = 0.5 * (r_grid[r_index] + r_grid[getIndex(r_i + 1, r_j, r_grid_width)]);
            int index_L = getIndex(i_center - 1, j_center, grid_width);
            grid[index_L] = r_grid[r_index];
            /*** fill top right corners and bottom left corners ***/
            if (r_j < (r_grid_height - 1)) {
                int index_RT = getIndex(i_center + 1, j_center + 1, grid_width);
                grid[index_RT] = 0.25 * (r_grid[r_index] + r_grid[getIndex(r_i + 1, r_j, r_grid_width)] + r_grid[getIndex(r_i, r_j + 1, r_grid_width)] + r_grid[getIndex(r_i + 1, r_j + 1, r_grid_width)]);
            }
            if (r_j > 0) {
                int index_BL = getIndex(i_center - 1, j_center - 1, grid_width);
                grid[index_BL] = 0.5 * (r_grid[r_index] + r_grid[getIndex(r_i, r_j - 1, r_grid_width)]);
            }
        }
        /** fill along right side blocks **/
        else if (r_i == (r_grid_width - 1)) {
            /*** fill right center side blocks, left center side blocks are already filled ***/
            int index_R = getIndex(i_center + 1, j_center, grid_width);
            grid[index_R] = r_grid[r_index];
            /*** fill top right corners, (except one) bottom left corners are already filled ***/
            if (r_j < (r_grid_height - 1)) {
                int index_TR = getIndex(i_center + 1, j_center + 1, grid_width);
                grid[index_TR] = 0.5 * (r_grid[r_index] + r_grid[getIndex(r_i, r_j + 1, r_grid_width)]);
            }
        }
        /** fill along to the bottom blocks **/
        else if (r_j == 0) {
            /*** fill top center side blocks and bottom center side blocks***/
            int index_T = getIndex(i_center, j_center + 1, grid_width);
            grid[index_T] = 0.5 * (r_grid[r_index] + r_grid[getIndex(r_i, r_j + 1, r_grid_width)]);
            int index_B = getIndex(i_center, j_center - 1, grid_width);
            grid[index_B] = r_grid[r_index];
            /*** fill top right corners and bottom left corners ***/
            if (r_i < (r_grid_width - 1)) {
                int index_RT = getIndex(i_center + 1, j_center + 1, grid_width);
                grid[index_RT] = 0.25 * (r_grid[r_index] + r_grid[getIndex(r_i + 1, r_j, r_grid_width)] + r_grid[getIndex(r_i, r_j + 1, r_grid_width)] + r_grid[getIndex(r_i + 1, r_j + 1, r_grid_width)]);
            }
            if (r_i > 0) {
                int index_BL = getIndex(i_center - 1, j_center - 1, grid_width);
                grid[index_BL] = 0.5 * (r_grid[r_index] + r_grid[getIndex(r_i - 1, r_j, r_grid_width)]);
            }
        }
        /** fill along top blocks **/
        else if (r_j == (r_grid_height - 1)) {
            /*** top center side blocks, bottom center side blocks are alredy filled ***/
            int index_T = getIndex(i_center, j_center + 1, grid_width);
            grid[index_T] = r_grid[r_index];
            /*** fill top right corners, bottom left corners are already filled (except one) ***/
            if (r_i < (r_grid_width - 1)) {
                int index_TR = getIndex(i_center + 1, j_center + 1, grid_width);
                grid[index_TR] = 0.5 * (r_grid[r_index] + r_grid[getIndex(r_i + 1, r_j, r_grid_width)]);
            }
        }
        /** fill bottom left BL corner block ***/
        if (r_i == 0 && r_j == 0) {
            int index_BL = getIndex(i_center - 1, j_center - 1, grid_width);
            grid[index_BL] = r_grid[r_index];
        }
        /** fill top right TR corner block **/
        else if (r_i == (r_grid_width - 1) && r_j == (r_grid_height - 1)) {
            int index_TR = getIndex(i_center + 1, j_center + 1, grid_width);
            grid[index_TR] = r_grid[r_index];
        }
        /** fill bottom right BR corner block **/
        else if (r_i == (r_grid_width - 1) && r_j == 0) {
            int index_BR = getIndex(i_center + 1, j_center - 1, grid_width);
            grid[index_BR] = r_grid[r_index];
        }
        /** fill top left TL corner block**/
        else if (r_i == 0 && r_j == (r_grid_height - 1)) {
            int index_TL = getIndex(i_center - 1, j_center + 1, grid_width);
            grid[index_TL] = r_grid[r_index];
        }
    }
    
}
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
