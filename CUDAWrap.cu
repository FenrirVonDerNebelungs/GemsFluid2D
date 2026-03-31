#include "CUDAWrap.h"
#include "Advection.cuh"
CUDAWrap::CUDAWrap(
    int blocks_side_dim,
    int threads_side_dim,
    double in_delta_t,
    double in_delta_x,
    double in_nu,
    int in_max_jacobi_loops
) : numBlocks_side(blocks_side_dim),
    numThreads_side(threads_side_dim),
    grid_width(blocks_side_dim * threads_side_dim),
    grid_height(blocks_side_dim * threads_side_dim),
    delta_t(in_delta_t),
    delta_x(in_delta_x),
    nu(in_nu),
    max_jacobi_loops(in_max_jacobi_loops)
{
    numBlocks_for_1D = blocks_side_dim * blocks_side_dim;
	numThreads_for_1D = threads_side_dim * threads_side_dim;
}
CUDAWrap::~CUDAWrap() {
    ;
}
void CUDAWrap::apply_force(double* Ux[], double* Uy[], int frame_index, s_force& force) {
	dim3 numBlocks(numBlocks_side, numBlocks_side);
	dim3 numThreads(numThreads_side, numThreads_side);
    s_frame_index frame_i = getFrameIndex(frame_index);
    double center_i = static_cast<double>(force.i);
    double center_j = static_cast<double>(force.j);
    double2 center = make_double2(center_i, center_j);
    double2 Force = make_double2(force.Fx_c, force.Fy_c);
    applyForce_Core << <numBlocks, numThreads >> > (
        Ux[frame_i.out],
        Uy[frame_i.out],
        Ux[frame_i.in],
        Uy[frame_i.in],
        center,
        Force,
        delta_t,
        force.inv_Rsqrd,
        grid_width);
    cudaError_t cudaStatus = cudaDeviceSynchronize();
}
void CUDAWrap::advection(double* Ux[], double* Uy[], int frame_index) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    s_frame_index frame_i = getFrameIndex(frame_index);
	advection_Core<< <numBlocks, numThreads >> > (
        Ux[frame_i.out], 
        Uy[frame_i.out], 
        Ux[frame_i.in], 
        Ux[frame_i.in], 
        delta_t, 
        delta_x, 
        grid_width, 
        grid_height);
	cudaError_t cudaStatus = cudaDeviceSynchronize();
}
void CUDAWrap::advection_backtrace(double* relPos_i, double* relPos_j, /*const*/ double* Ux[], /*const*/ double* Uy[], int frame_index) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    s_frame_index frame_i = getFrameIndex(frame_index);
    advection_backtrace_Core << < numBlocks, numThreads >> > (
        relPos_i,
        relPos_j,
        Ux[frame_i.in],
        Uy[frame_i.in],
        delta_t,
        delta_x,
        grid_width,
        grid_height
        );
    cudaError_t cudaStatus = cudaDeviceSynchronize();
}

void CUDAWrap::viscous_diffusion(double* Ux[], double* Uy[], double* scratch, int frame_index) {
	cudaError_t cudaStatus = cudaSuccess;
    /* jacobi with \alpha = frac{\deltax^2}{\nu \deltat} and b=\vec{u} and \beta = 4+\alpha */
    static double alpha = delta_x * delta_x / (nu * delta_t);
    static double beta = 4 + alpha;
    static double rbeta = 1.0 / beta;
    copy_memory << <numBlocks_for_1D, numThreads_for_1D >> > (scratch, Ux[frame_index]);
    cudaStatus = cudaDeviceSynchronize();
    jacobi_loop(Ux, scratch, frame_index, alpha, rbeta, false);
    copy_memory << <numBlocks_for_1D, numThreads_for_1D >> > (scratch, Uy[frame_index]);
    cudaStatus = cudaDeviceSynchronize();
    jacobi_loop(Uy, scratch, frame_index, alpha, rbeta, false);
}
void CUDAWrap::compute_pressure(double* Ux[], double* Uy[], double* p[], double* scratch, int frame_index, int p_frame_index) {
	divergence(scratch, Ux[frame_index], Uy[frame_index]);
    /* jacobi with \alpha = -\deltax^2 and b = \frac{1}{\deltat} \nabla \cdot \vec{u} and \beta = 4 */
    static double alpha = -delta_x * delta_x;
    static double rbeta = 0.25;
    jacobi_loop(p, scratch, p_frame_index, alpha, rbeta, Ux[frame_index], Uy[frame_index]);
}
void CUDAWrap::subtract_pressure_gradient(double* Ux[], double* Uy[], double* p[], int frame_index, int p_frame_index) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    s_frame_index frame_i = getFrameIndex(frame_index);
	static double inv_2delta_x = 1.0 / (2.0*delta_x);
    subtractGradient << <numBlocks, numThreads >> > (
        Ux[frame_i.out],
        Uy[frame_i.out],
        Ux[frame_i.in],
        Uy[frame_i.in],
        p[frame_i.in],
        inv_2delta_x,
        grid_width,
        grid_height);
	cudaError_t cudaStatus = cudaDeviceSynchronize();
}
void CUDAWrap::jacobi_frame(
    double* frame_out, 
    const double* frame_in, 
    const double* b, 
    const double* Wx,
    const double* Wy,
    const double& alpha, 
    const double& rbeta) 
{
    cudaError_t cudaStatus = cudaSuccess;
	dim3 numBlocks(numBlocks_side, numBlocks_side);
	dim3 numThreads(numThreads_side, numThreads_side);
    jacobi << <numBlocks, numThreads >> > (frame_out, frame_in, b, alpha, rbeta, grid_width, grid_height);
    cudaStatus = cudaDeviceSynchronize();
    if (Wx != nullptr)
        jacobi_boundary_pressure << <numBlocks_side, numThreads_side >> > (frame_out, Wx, Wy, frame_in, b, delta_x, grid_width, grid_height);
    cudaStatus = cudaDeviceSynchronize();
}
void CUDAWrap::jacobi_loop(double* X[], double* b, int frame_index, double alpha, double rbeta, const double* Wx, const double* Wy) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    cudaError_t cudaStatus = cudaSuccess;
    int original_frame_index = frame_index;
    s_frame_index frame_i = getFrameIndex(frame_index);
    int num_jacobi_loops = 0;
    do {
		jacobi_frame(X[frame_i.out], X[frame_i.in], b, Wx, Wy, alpha, rbeta);
        swapFrameIndexes(frame_i);
        num_jacobi_loops++;
    } while (num_jacobi_loops <= max_jacobi_loops);
    fixFramePointers(X, frame_i, original_frame_index);/* set Ux so that frame_out_index will point to the Ux results of the jacobi loop*/
}

void CUDAWrap::divergence(double* div, const double* Ux, const double* Uy) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
	static double inv_2delta_x = 1.0 / (2.0 * delta_x);
    divergence_Core << <numBlocks, numThreads >> > (div, Ux, Uy, inv_2delta_x, grid_width, grid_height);
    cudaError_t cudaStatus = cudaDeviceSynchronize();
}
void CUDAWrap::runFrame(double* Ux[], double* Uy[], double* p[], double* scratch, int& frame_index, int& p_frame_index, s_force& force) {
	int frame_in_index = frame_index;
    int p_frame_in_index = p_frame_index;
	advection(Ux, Uy, frame_in_index);
	reverseFrameIndex(frame_in_index);
	viscous_diffusion(Ux, Uy, scratch, frame_in_index);
	reverseFrameIndex(frame_in_index);
    if(force.active) {
        apply_force(Ux, Uy, frame_in_index, force);
		reverseFrameIndex(frame_in_index);
	}
	compute_pressure(Ux, Uy, p, scratch, frame_in_index, p_frame_in_index);
	reverseFrameIndex(p_frame_in_index);
	subtract_pressure_gradient(Ux, Uy, p, frame_in_index, p_frame_in_index);
    reverseFrameIndex(frame_in_index);
	frame_index = frame_in_index;
	p_frame_index = p_frame_in_index;
}
void CUDAWrap::bilinearAprox_scaledFrame(double* Ux_scaled, double* Uy_scaled, const double* Ux, const double* Uy, int scale_factor) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    bilinearAprox_scaledFrame_Core << <numBlocks, numThreads >> > (
        Ux_scaled,
        Uy_scaled,
        Ux,
        Uy,
        grid_width,
        grid_height,
        scale_factor);
    cudaError_t cudaStatus = cudaDeviceSynchronize();
}

int CUDAWrap::runNV(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames) {
    unsigned int size = grid_width * grid_height;
    double* dev_Ux[] = { 0,0 };/*two frames */
    double* dev_Uy[] = { 0,0 };
    double* dev_p[] = { 0,0 };
    double* scratch = 0;

    cudaError_t cudaStatus = cudaSuccess;
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaSetDevice failed!");

    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&(dev_Ux[0]), size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&(dev_Ux[1]), size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&(dev_Uy[0]), size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&(dev_Uy[1]), size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&(dev_p[0]), size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&(dev_p[1]), size * sizeof(double));
    if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&scratch, size * sizeof(double));
    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaMalloc failed!");


    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(dev_Ux[0], Ux, size * sizeof(double), cudaMemcpyHostToDevice);

    if (cudaStatus == cudaSuccess) {
        int frames_run = 0;
        int frame_index = 0;
        int p_frame_index = 0;
        do {
			runFrame(dev_Ux, dev_Ux, dev_p, scratch, frame_index, p_frame_index, force);
			frames_run++;
        } while (frames_run <= sim_frames);
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess)
            fprintf(stderr, "Cuda kernel launches failed:%s\n", cudaGetErrorString(cudaStatus));
    }
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(Ux, dev_Ux[0], size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(Uy, dev_Uy[0], size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(pressure, dev_p[0], size * sizeof(double), cudaMemcpyDeviceToHost);

    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaMemcpy failed!");

    cudaFree(dev_Ux[0]);
    cudaFree(dev_Ux[1]);
    cudaFree(dev_Uy[0]);
	cudaFree(dev_Uy[1]);
    cudaFree(dev_p[0]);
	cudaFree(dev_p[1]);
	cudaFree(scratch);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "runCUDA failed");
        return 1;
    }
	return 0;
}

int CUDAWrap::runCUDA(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames) {
    return runNV(Ux, Uy, pressure, force, sim_frames);
}
