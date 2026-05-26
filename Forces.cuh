#include "CUDABase.cuh"

__global__ void applyForce_Core(
    double* Ux_out,
    double* Uy_out,
    double* dye_out,
    const double* Ux_in,
    const double* Uy_in,
    const double* dye_in,
    const double2 center/*in i,j*/,
    const double2 F_c,
    const double2 dye_center,
    double dye_delta_mag,
    double delta_t,
    double inv_rsqrd,
    double inv_dye_rsqrd,
    const int grid_width)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int center_index = getIndex(i, j, grid_width);
    /* frac{(x-x_p)^2 + (y-y_p)^2}{r}*/
    double x_dist = (double)i - center.x;
    double y_dist = (double)j - center.y;
    double arg = -(x_dist * x_dist + y_dist * y_dist) * inv_rsqrd;
    double exp_arg = exp(arg);

    double dye_x_dist = (double)i - dye_center.x;
    double dye_y_dist = (double)j - dye_center.y;
    double dye_arg = -(dye_x_dist * dye_x_dist + dye_y_dist * dye_y_dist) * inv_dye_rsqrd;
    double dye_exp_arg = exp(dye_arg);

    double force_x = F_c.x * exp_arg;
    double force_y = F_c.y * exp_arg;
    double delta_ux = delta_t * force_x;
    double delta_uy = delta_t * force_y;
    Ux_out[center_index] = Ux_in[center_index] + delta_ux;
    Uy_out[center_index] = Uy_in[center_index] + delta_uy;

    double delta_dye = dye_delta_mag * dye_exp_arg;
    dye_out[center_index] = dye_in[center_index] + delta_dye;
}

__device__ void viscous_diffusion_for_component(
    double* frame_out,
    int i,
    int j,
    const double* frame_in,
    double nu/*viscosity constant*/,
    double inv_delta_Xsqrd,
    double delta_t,
    double grid_width,
    double grid_height)
{
    double xL = 0.0, xR = 0.0;
    double xT = 0.0, xB = 0.0;

    if (i > 0)
        xL = frame_in[getIndex(i - 1, j, grid_width)];
    if (i < (grid_width - 1))
        xR = frame_in[getIndex(i + 1, j, grid_width)];
    if (j > 0)
        xB = frame_in[getIndex(i, j - 1, grid_width)];
    if (j < (grid_height - 1))
        xT = frame_in[getIndex(i, j + 1, grid_width)];

    double center_val = frame_in[getIndex(i, j, grid_width)];
    double laplacian = (xL + xR + xT + xB - 4.0 * center_val) * inv_delta_Xsqrd;
    double new_velocity_component = center_val + nu * laplacian * delta_t;
    frame_out[getIndex(i, j, grid_width)] = new_velocity_component;
}
__global__ void viscous_diffusion_core(
    double* Ux_out,
    double* Uy_out,
    const double* Ux_in,
    const double* Uy_in,
    double nu/*viscosity constant*/,
    double inv_delta_Xsqrd,
    double delta_t,
    double grid_width,
    double grid_height)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    viscous_diffusion_for_component(Ux_out, i, j, Ux_in, nu, inv_delta_Xsqrd, delta_t, grid_width, grid_height);
    viscous_diffusion_for_component(Uy_out, i, j, Uy_in, nu, inv_delta_Xsqrd, delta_t, grid_width, grid_height);
}
__global__ void sliding_wall_core(
    double* Ux_in_out,
    double Ux)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    Ux_in_out[i] = Ux;
}