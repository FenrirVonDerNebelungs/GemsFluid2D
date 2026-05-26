#include "Pressure.cuh"

__device__ __constant__ int g_expansion_filter_half_wh=BASE_JACOBI_EXPANSION_FILTER_HALF_WH;


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
__global__ void Xgrid_reduction_dup(
    double* r_grid1,
    double* r_grid2,
    const int r_grid_width,
    const int r_grid_height,
    const double* grid1,
    const int grid_width,
    const int grid_height,
    const double* filter
)
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
    r_grid2[r_index] = r_grid1[r_index];
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

__device__ bool Xgrid_expansion_getCellValue(
	double& expanded_grid_cell_value,
    const double4& w_,
    const double2& I_,
    const double2& J_,
	const double* r_grid,
    const int r_grid_width,
    const int r_grid_height)
{
    bool retval = false;
    double4 corners = make_double4(0.0, 0.0, 0.0, 0.0);
    if (I_.x >= 0 && I_.y < r_grid_width && J_.x >= 0 && J_.y < r_grid_height) {
        /* L=x, R=y  B=x T=y*/
        /* BL=x, BR=y, TR=z, TL=w*/
        corners.x = r_grid[getIndex(I_.x, J_.x, r_grid_width)]; //BL
        corners.y = r_grid[getIndex(I_.y, J_.x, r_grid_width)]; //BR
        corners.z = r_grid[getIndex(I_.y, J_.y, r_grid_width)]; //TR
        corners.w = r_grid[getIndex(I_.x, J_.y, r_grid_width)]; //TL
        retval = true;
    }else
		retval = false;
    expanded_grid_cell_value = w_.x * corners.x + w_.y * corners.y + w_.z * corners.z + w_.w * corners.w;
    return retval;
}
__device__ void Xgrid_expansion_getCellValue_Edge(
    double& expanded_grid_cell_value,
    const double4& w_in,
    double l,
    double s,
    const double2& I_,
    const double2& J_,
    const double* r_grid,
    const int r_grid_width,
    const int r_grid_height)
{
	double4 w_ = w_in;
    double4 corners = make_double4(0.0, 0.0, 0.0, 0.0);
        /* L=x, R=y  B=x T=y*/
        /* BL=x, BR=y, TR=z, TL=w*/
    if (I_.x < 0) { //L
        if (J_.x < 0) //B
            w_.z = 1.0;//TR
        else if (J_.y >= r_grid_height) //T 
            w_.y = 1.0; //BR
        else {
            w_.z /= l; //TR
            w_.y /= l; //BR
        }
    }
    else if (I_.y >= r_grid_width) {
        if (J_.x < 0) //B
            w_.w = 1.0; //TL
        else if (J_.y >= r_grid_height) //T
            w_.x = 1.0; //BL
        else {
            w_.w /= l; //TL
            w_.x /= l; //BL
        }
    }
    else {/*r_I_L and r_I_R are ok */
        if (J_.x < 0) { //B
            w_.z /= l; //TR
            w_.w /= l; //TL
        }
        else { /*r_J_T must be >= r_grid_height*/
            w_.y /= l; //BR
            w_.x /= l; //BL
        }
    }

    if(I_.x>=0 && J_.x >=0)
        corners.x = r_grid[getIndex(I_.x, J_.x, r_grid_width)]; //BL
	if (I_.y < r_grid_width && J_.x >= 0)
        corners.y = r_grid[getIndex(I_.y, J_.x, r_grid_width)]; //BR
	if (I_.y < r_grid_width && J_.y < r_grid_height)
        corners.z = r_grid[getIndex(I_.y, J_.y, r_grid_width)]; //TR
	if (I_.x >= 0 && J_.y < r_grid_height)
        corners.w = r_grid[getIndex(I_.x, J_.y, r_grid_width)]; //TL
    
    expanded_grid_cell_value = w_.x * corners.x + w_.y * corners.y + w_.z * corners.z + w_.w * corners.w;
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
	double4 w_ = make_double4(0.0, 0.0, 0.0, 0.0);
	double2 I_ = make_double2(0.0, 0.0);
	double2 J_ = make_double2(0.0, 0.0);
    double4 corners = make_double4(0.0, 0.0, 0.0, 0.0);
    double expanded_grid_cell_value = 0.0;
    /* BL=x, BR=y, TR=z, TL=w*/
    if (i % 2 == 0) {
        I_.y = i / 2; //R
        w_.z = l; //TR
        w_.y = l; //BR
        I_.x = I_.y - 1; //L
        w_.w = s; //TL
        w_.x = s; //BL
    }
    else {
        I_.x = i / 2; //L
		w_.w = l; //TL
        w_.x = l; //BL
        I_.y = I_.x + 1; //R
        w_.z = s; //TR
        w_.y = s; //BR
    }
    if (j % 2 == 0) {
        J_.y = j / 2; //T
        w_.w *= l; //TL
        w_.z *= l; //TR
        J_.x = J_.y - 1; //B
        w_.x *= s; //BL
        w_.y *= s; //BR
    }
    else {
        J_.x = j / 2; //B
        w_.x *= l; //BL
        w_.y *= l; //BR
        J_.y = J_.x + 1; //T
        w_.w *= s; //TL
        w_.z *= s; //TR
    }
    if( !Xgrid_expansion_getCellValue(expanded_grid_cell_value,w_, I_, J_, r_grid, r_grid_width, r_grid_height) )
		Xgrid_expansion_getCellValue_Edge(expanded_grid_cell_value, w_, l, s, I_, J_, r_grid, r_grid_width, r_grid_height);
    int grid_index = getIndex(i, j, grid_width);
    grid[grid_index] = expanded_grid_cell_value;
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
