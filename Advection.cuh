#include "Jacobi.cuh"
__device__ double2 getRelBackTracedPosition(
    double2 U,
    const double delta_t,
    const double delta_x) 
{
    double Ux_back = -U.x * delta_t / delta_x;
    double Uy_back = -U.y * delta_t / delta_x;
    double2 relPos = make_double2(Ux_back, Uy_back);
    return relPos;
}
__device__ double2 getCorrectedOutOfRangeRelBackTracedPosition(
    const double2 relPos,
    const int i,
    const int j,
    const int grid_width,
    const int grid_height) 
{
    double2 newRelPos = make_double2(relPos.x, relPos.y);
    double distToXwall = (relPos.x < 0) ? i : (grid_width - 1) - i;
    double distToYwall = (relPos.y < 0) ? j : (grid_height - 1) - j;
    double abs_x_relPos = abs(relPos.x);
    double abs_y_relPos = abs(relPos.y);
    double2 x_edge_cross_point = make_double2(0.0, 0.0);
    if (distToXwall < abs_x_relPos && relPos.x != 0.0) {
        /* went over an x wall and abs_x_relPos>0*/
        double slope = relPos.y / relPos.x;
        x_edge_cross_point.x = (relPos.x < 0) ? -distToXwall : distToXwall;
        x_edge_cross_point.y = slope * x_edge_cross_point.x;
        newRelPos = x_edge_cross_point;
    }
    double2 y_edge_cross_point = make_double2(0.0, 0.0);
    if (distToYwall < abs_y_relPos && relPos.y != 0.0) {
        double slope = relPos.x / relPos.y;
        y_edge_cross_point.y = (relPos.y < 0) ? -distToYwall : distToYwall;
        y_edge_cross_point.x = slope * y_edge_cross_point.y;
        newRelPos = y_edge_cross_point;
    }
    if (distToXwall < abs_x_relPos && distToYwall < abs_y_relPos) {
        if (relPos.x != 0.0 && relPos.y != 0.0) {
            double dist_x_edge_cross_point_to_y_edge = distToYwall - abs(x_edge_cross_point.y);
            double dist_y_edge_cross_point_to_x_edge = distToXwall - abs(y_edge_cross_point.x);
            if (dist_x_edge_cross_point_to_y_edge > dist_y_edge_cross_point_to_x_edge)
                newRelPos = x_edge_cross_point; /* else newRelPos is already set to y_edge_cross_point*/
        } /* if relPos.x!=0 relPos.y==0 abs_y_relPos==0 so distToYwall>= abs_y_relPos
           * if relPos.y!=0 .. likewise... */
    }
    return newRelPos;
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
    s_four_corners relCorners = getRelFourCorners(i, j, relPos);
    s_four_corners corners = addCornersTo_ij(i, j, relCorners);
    double2 U_out;
    if (inRange_corners(corners, grid_width, grid_height)) {
        /* find the interpolated values for U_x, U_y*/
        U_out = bilinearAprox_Core(Ux_in, Uy_in, relPos, relCorners, corners, delta_x, grid_width);
    }
    else {
        double2 shortened_relPos = getCorrectedOutOfRangeRelBackTracedPosition(relPos, i, j, grid_width, grid_height);
        s_four_corners edge_relCorners = getRelFourCorners(i, j, shortened_relPos);
        s_four_corners edge_corners = addCornersTo_ij(i, j, edge_relCorners);
        U_out = bilinearAprox_Core(Ux_in, Uy_in, shortened_relPos, edge_relCorners, edge_corners, delta_x, grid_width);
    }

    Ux_out[center_index] = U_out.x;
    Uy_out[center_index] = U_out.y;
}