#include "CUDAWrap.h"
#include "Advection.cuh"
CUDAWrap::CUDAWrap(
    int blocks_side_dim,
    int threads_side_dim,
    double in_delta_t,
    double in_delta_x,
    double in_nu,
	int jacobi_minBlocks_side_dim,
	int jacobi_minThreads_side_dim,
    int in_max_jacobi_loops
) : numBlocks_side(blocks_side_dim),
    numThreads_side(threads_side_dim),
    grid_width(blocks_side_dim * threads_side_dim),
    grid_height(blocks_side_dim * threads_side_dim),
    delta_t(in_delta_t),
    delta_x(in_delta_x),
    nu(in_nu),
	jacobi_minBlocks_side(jacobi_minBlocks_side_dim),
	jacobi_minThreads_side(jacobi_minThreads_side_dim),
    max_jacobi_loops(in_max_jacobi_loops)
{
    numBlocks_for_1D = blocks_side_dim * blocks_side_dim;
	numThreads_for_1D = threads_side_dim * threads_side_dim;
    if(numBlocks_side% jacobi_minBlocks_side != 0 || numThreads_side % jacobi_minThreads_side != 0)
		fprintf(stderr, "numBlocks and numThreads must be a factor of 2 times jacobi_minBlocks_side and jacobi_minThreads_side respectively ");
    int block_jacobi_expansion_factor = numBlocks_side / jacobi_minBlocks_side;
	int thread_jacobi_expansion_factor = numThreads_side / jacobi_minThreads_side;
	if (block_jacobi_expansion_factor % 2 != 0 || thread_jacobi_expansion_factor % 2 != 0)
		fprintf(stderr, "numBlocks and numThreads must be a factor of 2 times jacobi_minBlocks_side and jacobi_minThreads_side respectively ");
	jacobi_block_expansion_pow = findPow2(block_jacobi_expansion_factor);
	jacobi_thread_expansion_pow = findPow2(thread_jacobi_expansion_factor);
	jacobi_stack_height = jacobi_block_expansion_pow + jacobi_thread_expansion_pow;
	jacobi_scratch_stack_sizes = nullptr;
    jacobi_scratch_stack = nullptr;
	jacobi_scratch = nullptr;
	b_stack = nullptr;
    if (jacobi_stack_height > 0) {
		jacobi_scratch_stack_sizes = new int[jacobi_stack_height];
        jacobi_scratch_stack = new double* [jacobi_stack_height];
		b_stack = new double* [jacobi_stack_height];
        int mul_2 = 1;
        for(int i=0; i<jacobi_stack_height; i++) {
			jacobi_scratch_stack_sizes[i] = mul_2 * jacobi_minBlocks_side * jacobi_minThreads_side;
			mul_2 *= 2;
			jacobi_scratch_stack[i] = nullptr;
		}
    }
}
CUDAWrap::~CUDAWrap() {
    if (jacobi_scratch_stack != nullptr) {
		delete[] jacobi_scratch_stack;
    }
    if (jacobi_scratch_stack_sizes != nullptr) {
        delete[] jacobi_scratch_stack_sizes;
    }
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
        Uy[frame_i.in], 
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
void CUDAWrap::jacobi_b_stack(const double* b) {
    int jacobi_stack_max_index = jacobi_stack_height - 1;
    copy_memory << <numBlocks_for_1D, numThreads_for_1D >> > (b_stack[jacobi_stack_max_index], b);
    int r_grid_width = grid_width;
    int r_grid_height = grid_height;
    int base_grid_width = grid_width;
    int base_grid_height = grid_height;
    int red_factor = 1;
    for (int i = 0; i < jacobi_stack_max_index; i++)
    {
        red_factor *= 2;
        r_grid_width /= 2;
        r_grid_height /= 2;
        int stack_base_i = jacobi_stack_max_index - i;
        int r_stack_i = jacobi_stack_max_index - i - 1;
        int numBlocks_s=0, int numThreads_s=0;
        if (!find_reduced_BlocksNThreads(numBlocks_s, numThreads_s, red_factor))
            break;
        dim3 numBlocks(numBlocks_s, numBlocks_s);
        dim3 numThreads(numThreads_s, numThreads_s);
        Xgrid_reduction << <numBlocks, numThreads >> > (
            b_stack[r_stack_i],
            r_grid_width,
            r_grid_height,
            b_stack[stack_base_i],
            base_grid_width,
            base_grid_height);
        cudaError_t cudaStatus = cudaDeviceSynchronize();
        base_grid_width /= 2;
        base_grid_height /= 2;
    }
}
void CUDAWrap::jacobi_run_stack(
    const double* frame_in,
    const double& alpha,
    const double& rbeta)
{
	int jacobi_stack_max_index = jacobi_stack_height - 1;
    copy_memory << <numBlocks_for_1D, numThreads_for_1D >> > (jacobi_scratch_stack[jacobi_stack_max_index], frame_in);
    int r_grid_width = grid_width;
	int r_grid_height = grid_height;
	int base_grid_width = grid_width;
	int base_grid_height = grid_height;
    int red_factor = 1;
    for (int i = 0; i < jacobi_stack_max_index; i++)
    {
		red_factor *= 2;
        r_grid_width /= 2;
		r_grid_height /= 2;
		int stack_base_i = jacobi_stack_max_index - i;
        int r_stack_i = jacobi_stack_max_index - i - 1;
        int numBlocks_s=0, int numThreads_s=0;
        if (!find_reduced_BlocksNThreads(numBlocks_s, numThreads_s, red_factor))
            break;
        dim3 numBlocks(numBlocks_s, numBlocks_s);
        dim3 numThreads(numThreads_s, numThreads_s);
		Xgrid_reduction << <numBlocks, numThreads >> > (
            jacobi_scratch_stack[r_stack_i], 
            r_grid_width, 
            r_grid_height, 
            jacobi_scratch_stack[stack_base_i],
            base_grid_width, 
            base_grid_height);
		cudaError_t cudaStatus = cudaDeviceSynchronize();
		base_grid_width /= 2;
		base_grid_height /= 2;
    }
    for (int r_stack_i = 0; r_stack_i < jacobi_stack_max_index; r_stack_i++) {
        int numBlocks_s=0, int numThreads_s=0;
        if (!find_reduced_BlocksNThreads(numBlocks_s, numThreads_s, red_factor))
            break;
        dim3 numBlocks(numBlocks_s, numBlocks_s);
        dim3 numThreads(numThreads_s, numThreads_s);
        jacobi << <numBlocks, numThreads >> > (
            jacobi_scratch,
            jacobi_scratch_stack[r_stack_i],
            b_stack[r_stack_i],
            alpha,
            rbeta,
            r_grid_width,
			r_grid_height);
        cudaError_t cudaStatus = cudaDeviceSynchronize();
        Xgrid_expansion << <numBlocks, numThreads >> > (
            jacobi_scratch_stack[r_stack_i + 1],
            r_grid_width * 2,
            r_grid_height * 2,
            jacobi_scratch,
            r_grid_width,
			r_grid_height);
		cudaStatus = cudaDeviceSynchronize();
		r_grid_width *= 2;
		r_grid_height *= 2;
		red_factor /= 2;
    }
}
void CUDAWrap::jacobi_base(
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
        jacobi_boundary_pressure << <numBlocks_side, numThreads_side >> > (frame_out, Wx, Wy, frame_in, b, alpha, rbeta, delta_x, grid_width, grid_height);
    cudaStatus = cudaDeviceSynchronize();
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
    if (jacobi_stack_height > 0) {
		jacobi_b_stack(b);
        jacobi_run_stack(frame_in, alpha, rbeta);
		int jacobi_stack_max_index = jacobi_stack_height - 1;
        jacobi_base(frame_out,jacobi_scratch_stack[jacobi_stack_max_index], b, Wx, Wy, alpha, rbeta);
    }else
		jacobi_base(frame_out, frame_in, b, Wx, Wy, alpha, rbeta);
}
void CUDAWrap::jacobi_loop(double* X[], double* b, int frame_index, double alpha, double rbeta, const double* Wx, const double* Wy) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    //cudaError_t cudaStatus = cudaSuccess;
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
    double* dev_Ux_0 = 0;
	double* dev_Ux_1 = 0;
    double* dev_Uy[] = { 0,0 };
	double* dev_Uy_0 = 0;
	double* dev_Uy_1 = 0;
    double* dev_p[] = { 0,0 };
	double* dev_p_0 = 0;
	double* dev_p_1 = 0;
    double* scratch = 0;

    cudaError_t cudaStatus = cudaSuccess;
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaSetDevice failed!");

    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&dev_Ux_0, size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&dev_Ux_1, size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&dev_Uy_0, size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&dev_Uy_1, size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&dev_p_0, size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&dev_p_1, size * sizeof(double));
    if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&scratch, size * sizeof(double));

    for (int i = 0; i < jacobi_stack_height; i++) {
        jacobi_scratch_stack[i] = 0;
		b_stack[i] = 0;
		if (cudaStatus == cudaSuccess)
			cudaStatus = cudaMalloc((void**)&jacobi_scratch_stack[i], jacobi_scratch_stack_sizes[i] * sizeof(double));
		if (cudaStatus == cudaSuccess)
			cudaStatus = cudaMalloc((void**)&b_stack[i], jacobi_scratch_stack_sizes[i] * sizeof(double));
    }
    jacobi_scratch = 0;
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMalloc((void**)&jacobi_scratch, size * sizeof(double));

    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaMalloc failed!");


    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(dev_Ux_0, Ux, size * sizeof(double), cudaMemcpyHostToDevice);
    if(cudaStatus==cudaSuccess)
		cudaStatus = cudaMemcpy(dev_Uy_0, Uy, size * sizeof(double), cudaMemcpyHostToDevice);

	dev_Ux[0] = dev_Ux_0;
	dev_Ux[1] = dev_Ux_1;
	dev_Uy[0] = dev_Uy_0;
	dev_Uy[1] = dev_Uy_1;
	dev_p[0] = dev_p_0;
	dev_p[1] = dev_p_1;

    if (cudaStatus == cudaSuccess) {
        int frames_run = 0;
        int frame_index = 0;
        int p_frame_index = 0;
        do {
			runFrame(dev_Ux, dev_Uy, dev_p, scratch, frame_index, p_frame_index, force);
			frames_run++;
        } while (frames_run <= sim_frames);
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess)
            fprintf(stderr, "Cuda kernel launches failed:%s\n", cudaGetErrorString(cudaStatus));
    }
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(Ux, dev_Ux_0, size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(Uy, dev_Uy_0, size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(pressure, dev_p_0, size * sizeof(double), cudaMemcpyDeviceToHost);

    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaMemcpy failed!");

	cudaFree(jacobi_scratch);
	for (int i = 0; i < jacobi_stack_height; i++) {
        cudaFree(jacobi_scratch_stack[i]);
        cudaFree(b_stack[i]);
	}
    cudaFree(dev_Ux_0);
    cudaFree(dev_Ux_1);
    cudaFree(dev_Uy_0);
	cudaFree(dev_Uy_1);
    cudaFree(dev_p_0);
	cudaFree(dev_p_1);
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

int CUDAWrap::findPow2(int val) {
    if (val <= 0)
        return 0;
    int pow = 0;
    while(val%2==0)
    {
        pow++;
		val /= 2;
    }
    return pow;
}