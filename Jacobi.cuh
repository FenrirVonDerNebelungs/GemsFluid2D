#include "Pressure.cuh"

__device__ __constant__ int g_expansion_filter_half_wh=BASE_JACOBI_EXPANSION_FILTER_HALF_WH;

__device__ int getReducedIndex(const int i, const int j, const int r_grid_width) {
	int reduction_i = i / 2;
	int reduction_j = j / 2;
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
__device__ double Xgrid_reduction_filter(
    const double* grid,
    const int grid_width,
    const int grid_height,
    const int i_lo_corner,
    const int j_lo_corner,
    const double* filter)
{
    double red_grid_value = 0.0;
    const int filter_wh = 2 * g_expansion_filter_half_wh;
    for (int j_filter = 0; j_filter < filter_wh; j_filter++) {
        int j = j_lo_corner + j_filter;
        for (int i_filter = 0; i_filter < filter_wh; i_filter++) {
            int i = i_lo_corner + i_filter;
            int grid_index = getIndex(i, j, grid_width);
            double grid_value = grid[grid_index];
            int filter_index = getIndex(i_filter, j_filter, filter_wh);
            double filter_weight = filter[filter_index];
            red_grid_value += filter_weight * grid_value;
        }
    }
    return red_grid_value;
}
__device__ double Xgrid_reduction_filter_edge(
    const double* grid,
    const int grid_width,
    const int grid_height,
    const int i_lo_corner,
    const int j_lo_corner,
    const double* filter)
{
    double red_grid_value = 0.0;
    double total_weight = 0.0;
    const int filter_wh = 2 * g_expansion_filter_half_wh;
    for (int j_filter = 0; j_filter < filter_wh; j_filter++) {
        int j = j_lo_corner + j_filter;
        for (int i_filter = 0; i_filter < filter_wh; i_filter++) {
            int i = i_lo_corner + i_filter;
            int grid_index = getGoodIndex(i, j, grid_width, grid_height);
            if (grid_index >= 0) {
                double grid_value = grid[grid_index];
                int filter_index = getIndex(i_filter, j_filter, filter_wh);
                double filter_weight = filter[filter_index];
                red_grid_value += filter_weight * grid_value;
                total_weight += filter_weight;
            }
        }
    }
    if (total_weight > 0.0)
        red_grid_value /= total_weight;
    return red_grid_value;
}

__device__ void Xgrid_reduction(
    double* r_grid,
    const int r_index,
    const double* grid,
    const int grid_width,
    const int grid_height,
    const int i_lo_corner,
    const int j_lo_corner,
    const double* filter) 
{
    const int filter_corner_steps_forward = 2*g_expansion_filter_half_wh-1;
    int i_hi_corner = i_lo_corner + filter_corner_steps_forward;
    int j_hi_corner = j_lo_corner + filter_corner_steps_forward;
    double red_grid_value = 0.0;
    if (i_lo_corner >= 0 && j_lo_corner >= 0 && i_hi_corner < grid_width && j_hi_corner < grid_height)
        red_grid_value = Xgrid_reduction_filter(grid, grid_width, grid_height, i_lo_corner, j_lo_corner, filter);
    else
        red_grid_value = Xgrid_reduction_filter_edge(grid, grid_width, grid_height, i_lo_corner, j_lo_corner, filter);
    r_grid[r_index] = red_grid_value;
}
__global__ void Xgrid_reduction_x4(
    double* r_grid1,
    double* r_grid2,
    double* r_grid3,
    double* r_grid4,
    const int r_grid_width,
    const int r_grid_height,
    const double* grid1,
    const double* grid2,
    const double* grid3,
    const double* grid4,
    const int grid_width,
    const int grid_height,
    const double* filter)
{
    const int filter_corner_steps_back = g_expansion_filter_half_wh - 1;
    int r_i = blockIdx.x * blockDim.x + threadIdx.x;
    int r_j = blockIdx.y * blockDim.y + threadIdx.y;
    int r_index = getIndex(r_i, r_j, r_grid_width);
    int i = 2 * r_i;
    int j = 2 * r_j;
    int i_lo_corner = i - filter_corner_steps_back;
    int j_lo_corner = j - filter_corner_steps_back;

    Xgrid_reduction(r_grid1, r_index, grid1, grid_width, grid_height, i_lo_corner, j_lo_corner, filter);
    Xgrid_reduction(r_grid2, r_index, grid2, grid_width, grid_height, i_lo_corner, j_lo_corner, filter);
    Xgrid_reduction(r_grid3, r_index, grid3, grid_width, grid_height, i_lo_corner, j_lo_corner, filter);
    Xgrid_reduction(r_grid4, r_index, grid4, grid_width, grid_height, i_lo_corner, j_lo_corner, filter);
}
/*r_grid is 1/2 width(height) of grid
* i=2I+1, j=2J+1, for grid cell at center of each reduced grid cell
*/
__device__ void Xgrid_reduction_paper(
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
__global__ void Xgrid_reduction_x4_paper(
    double* r_grid1,
	double* r_grid2,
	double* r_grid3,
	double* r_grid4,
    const int r_grid_width, 
    const int r_grid_height, 
    const double* grid1,
    const double* grid2,
	const double* grid3,
	const double* grid4,
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
	Xgrid_reduction_paper(r_grid1, r_index, grid1, i_center, i_L, i_R, i_B, i_T);
	Xgrid_reduction_paper(r_grid2, r_index, grid2, i_center, i_L, i_R, i_B, i_T);
	Xgrid_reduction_paper(r_grid3, r_index, grid3, i_center, i_L, i_R, i_B, i_T);
	Xgrid_reduction_paper(r_grid4, r_index, grid4, i_center, i_L, i_R, i_B, i_T);
}
/*grid must be 2width and 2height of r_grid
* blocks/threads are run over expanded grid*/
__global__ void Xgrid_expansion(
    double* grid,
    const int grid_width,
    const int grid_height,
    const double* r_grid,
    const int r_grid_width,
    const int r_grid_height)
{
    const double s = 0.25;
    const double l = 0.75;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int r_I_L = 0, r_I_R=0;
    int r_J_B=0, r_J_T=0;
    double w_TL = 0.0, w_TR = 0.0;
    double w_BL = 0.0, w_BR = 0.0;
    double val_TL = 0.0, val_TR = 0.0;
    double val_BL = 0.0, val_BR = 0.0;
    if (i % 2 == 0) {
        r_I_R = i / 2;
        w_TR = l;
        w_BR = l;
        r_I_L = r_I_R - 1;
        w_TL = s;
        w_BL = s;
    }
    else {
        r_I_L = i / 2;
        w_TL = l;
        w_BL = l;
        r_I_R = r_I_L + 1;
        w_TR = s;
        w_BR = s;
    }
    if (j % 2 == 0) {
        r_J_T = j / 2;
        w_TL *= l;
        w_TR *= l;
        r_J_B = r_J_T - 1;
        w_BL *= s;
        w_BR *= s;
    }
    else {
        r_J_B = j / 2;
        w_BL *= l;
        w_BR *= l;
        r_J_T = r_J_B + 1;
        w_TL *= s;
        w_TR *= s;
    }
    if (r_I_L >= 0 && r_I_R < r_grid_width && r_J_B >= 0 && r_J_T < r_grid_height) {
        val_BL = r_grid[getIndex(r_I_L, r_J_B, r_grid_width)];
        val_BR = r_grid[getIndex(r_I_R, r_J_B, r_grid_width)];
        val_TR = r_grid[getIndex(r_I_R, r_J_T, r_grid_width)];
        val_TL = r_grid[getIndex(r_I_L, r_J_T, r_grid_width)];
    }
    else {
        if (r_I_L < 0) {
            if (r_J_B < 0) {
                val_TR = r_grid[getIndex(r_I_R, r_J_T, r_grid_width)];
                w_TR = 1.0;
            }
            else if (r_J_T >= r_grid_height) {
                val_BR = r_grid[getIndex(r_I_R, r_J_B, r_grid_width)];
                w_BR = 1.0;
            }
            else {
                w_TR /= l;
                w_BR /= l;
                val_TR = r_grid[getIndex(r_I_R, r_J_T, r_grid_width)];
                val_BR = r_grid[getIndex(r_I_R, r_J_B, r_grid_width)];
            }
        }
        else if (r_I_R >= r_grid_width) {
            if (r_J_B < 0) {
                val_TL = r_grid[getIndex(r_I_L, r_J_T, r_grid_width)];
                w_TL = 1.0;
            }
            else if (r_J_T >= r_grid_height) {
                val_BL = r_grid[getIndex(r_I_L, r_J_B, r_grid_width)];
                w_BL = 1.0;
            }
            else {
                w_TL /= l;
                w_BL /= l;
                val_TL = r_grid[getIndex(r_I_L, r_J_T, r_grid_width)];
                val_BL = r_grid[getIndex(r_I_L, r_J_B, r_grid_width)];
            }
        }
        else {/*r_I_L and r_I_R are ok */
            if (r_J_B < 0) {
                w_TR /= l;
                w_TL /= l;
                val_TR = r_grid[getIndex(r_I_R, r_J_T, r_grid_width)];
                val_TL = r_grid[getIndex(r_I_L, r_J_T, r_grid_width)];
            }
            else { /*r_J_T must be >= r_grid_height*/
                w_BR /= l;
                w_BL /= l;
                val_BR = r_grid[getIndex(r_I_R, r_J_B, r_grid_width)];
                val_BL = r_grid[getIndex(r_I_L, r_J_B, r_grid_width)];
            }
        }
    }
    double expanded_grid_cell_value = w_BL * val_BL + w_BR * val_BR + w_TR * val_TR + w_TL * val_TL;
    int grid_index = getIndex(i, j, grid_width);
    grid[grid_index] = expanded_grid_cell_value;
}
__global__ void Xgrid_expansion_paper(
    double* grid, 
    const int grid_width, 
    const int grid_height, 
    const double* r_grid, 
    const int r_grid_width, 
    const int r_grid_height)
{ 
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    double div_cnt_I = 1.0;
    double div_cnt_J = 1.0;
    /* i=2I+1
     * j=2J+1 */
	int _I = (i - 1) / 2;
	int _J = (j - 1) / 2;
    int r_I = _I;
    int r_J = _J;
    int r_I_lo = -1, r_I_hi = -1;
	int r_J_lo = -1, r_J_hi = -1;
    if ((i - 1) % 2 == 0) {
        r_I = _I;
    }
    else {
        if (i > 0 && i < (grid_width - 1)) {
            r_I_lo = _I;
			r_I_hi = _I + 1;
            div_cnt_I+=1.0;
        }
        else {
            r_I = _I;
        }
    }
    if ((j - 1) % 2 == 0) {
        r_J = _J;
    }
    else {
        if (j > 0 && j < (grid_height - 1)) {
			r_J_lo = _J;
            r_J_hi = _J + 1;
            div_cnt_J+=1.0;
        }
        else {
            r_J = _J;
        }
    }
    double div_cnt = div_cnt_I * div_cnt_J;
    double val_sum = 0.0;
    if (div_cnt < 1.1) {
		val_sum = r_grid[getIndex(r_I, r_J, r_grid_width)];
    }
    else if (div_cnt < 2.1) {
        if (r_I_lo >= 0) {
			val_sum += r_grid[getIndex(r_I_lo, r_J, r_grid_width)];
			val_sum += r_grid[getIndex(r_I_hi, r_J, r_grid_width)];
        }
        else {
			val_sum += r_grid[getIndex(r_I, r_J_lo, r_grid_width)]; 
			val_sum += r_grid[getIndex(r_I, r_J_hi, r_grid_width)];
        }
    }
    else {
        val_sum += r_grid[getIndex(r_I_lo, r_J_lo, r_grid_width)];
        val_sum += r_grid[getIndex(r_I_lo, r_J_hi, r_grid_width)];
        val_sum += r_grid[getIndex(r_I_hi, r_J_lo, r_grid_width)];
		val_sum += r_grid[getIndex(r_I_hi, r_J_hi, r_grid_width)];
    }
	grid[getIndex(i, j, grid_width)] = val_sum / div_cnt;
}
/*grid must be 2*width and 2*height of r_grid
* blocks/threads run over REDUCED grid 
* should give same results as above */
__global__ void Xgrid_expansion_from_reduced(double* grid, const int grid_width, const int grid_height, const double* r_grid, const int r_grid_width, const int r_grid_height)
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
    const double* Wx_in,
	const double* Wy_in,
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
	const double g_mult_constant = 2.0*delta_x;
    double g_mult_constant_x = 0.0;
    double g_mult_constant_y = 0.0;
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
		g_mult_constant_y = g_mult_constant;
        xB = xT;
    }
    if(jT<0){
        g_index = getIndex(i, j_max, grid_width);
        g_mult_constant_y = -g_mult_constant;
        xT = xB;
    }
    if(iL<0){
        g_index = getIndex(0, j, grid_width);
		g_mult_constant_x = g_mult_constant;
        xL = xR;
    }
    if(iR<0){
        g_index = getIndex(i_max, j, grid_width);
        g_mult_constant_x = -g_mult_constant;
        xR = xL;
    }
    g = g_mult_constant_x * Wx_in[g_index] + g_mult_constant_y*Wy_in[g_index];
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
    
    if (k < (grid_width - 1)) {
        j = k;
        jB = j - 1;
        jT = j + 1;
        /* go up i=0 */
        i = 0;
        iR = i + 1;
        iL = -1;
        jacobi_boundary_pressure_doSide(frame_out, frame_in, Wx_in, Wy_in, b_in, i, iL, iR, j, jT, jB, alpha, rbeta, delta_x, grid_width, grid_height);
        /* go up i=i_max-1 */
        i = grid_width - 1;
        iR = -1;
        iL = i - 1;
        jacobi_boundary_pressure_doSide(frame_out, frame_in, Wx_in, Wy_in, b_in, i, iL, iR, j, jT, jB, alpha, rbeta, delta_x, grid_width, grid_height);

        i = k;
        iL = i - 1;
        iR = i + 1;
        /* go along j=0 */
        j = 0;
        jT = j + 1;
        jB = -1;
        jacobi_boundary_pressure_doSide(frame_out, frame_in, Wx_in, Wy_in, b_in, i, iL, iR, j, jT, jB, alpha, rbeta, delta_x, grid_width, grid_height);
        /* go along j=j_max-1 */
        j = grid_height - 1;
        jT = -1;
        jB = j - 1;
        jacobi_boundary_pressure_doSide(frame_out, frame_in, Wx_in, Wy_in, b_in, i, iL, iR, j, jT, jB, alpha, rbeta, delta_x, grid_width, grid_height);
    }
    else {
        j = grid_height - 1;
        jB = j - 1;
        jT = -1;
        i = 0;
        iR = i + 1;
		iL = -1;
        jacobi_boundary_pressure_doSide(frame_out, frame_in, Wx_in, Wy_in, b_in, i, iL, iR, j, jT, jB, alpha, rbeta, delta_x, grid_width, grid_height);
		i = grid_width - 1;
        iR = -1;
        iL = i - 1;
		jacobi_boundary_pressure_doSide(frame_out, frame_in, Wx_in, Wy_in, b_in, i, iL, iR, j, jT, jB, alpha, rbeta, delta_x, grid_width, grid_height);
        j = 0;
		jT = j + 1;
        jB = -1;
		jacobi_boundary_pressure_doSide(frame_out, frame_in, Wx_in, Wy_in, b_in, i, iL, iR, j, jT, jB, alpha, rbeta, delta_x, grid_width, grid_height);
    }
}
