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
    int in_max_jacobi_loops,
	int in_max_jacobi_force_loops
) : numBlocks_side(blocks_side_dim),
    numThreads_side(threads_side_dim),
    grid_width(blocks_side_dim * threads_side_dim),
    grid_height(blocks_side_dim * threads_side_dim),
    delta_t(in_delta_t),
    delta_x(in_delta_x),
    nu(in_nu),
	jacobi_minBlocks_side(jacobi_minBlocks_side_dim),
	jacobi_minThreads_side(jacobi_minThreads_side_dim),
    m_max_jacobi_loops(in_max_jacobi_loops), 
	m_max_jacobi_force_loops(in_max_jacobi_force_loops),
	m_num_jacobi_loops(in_max_jacobi_loops)
{
    for(int i=0; i<2; i++){
        m_dev_Ux[i] = nullptr;
        m_dev_Uy[i] = nullptr;
		m_dev_p0[i] = nullptr;
        m_dev_p[i] = nullptr;
		m_dev_dye[i] = nullptr;
	}
	m_dev_scratch = nullptr;
    m_filter = nullptr;
    for(int i=0; i<4; i++) {
        advection_indexes[i] = nullptr;
        advection_dist[i] = nullptr;
	}

    numBlocks_for_1D = numBlocks_side * numBlocks_side;
	numThreads_for_1D = numThreads_side * numThreads_side;
    if(numBlocks_side% jacobi_minBlocks_side != 0 || numThreads_side % jacobi_minThreads_side != 0)
		fprintf(stderr, "numBlocks and numThreads must be a factor of 2 times jacobi_minBlocks_side and jacobi_minThreads_side respectively ");
    int block_jacobi_expansion_factor = numBlocks_side / jacobi_minBlocks_side;
	int thread_jacobi_expansion_factor = numThreads_side / jacobi_minThreads_side;
	if (block_jacobi_expansion_factor % 2 != 0 || thread_jacobi_expansion_factor % 2 != 0)
		fprintf(stderr, "numBlocks and numThreads must be a factor of 2 times jacobi_minBlocks_side and jacobi_minThreads_side respectively ");
	jacobi_block_expansion_pow = findLog_base2(block_jacobi_expansion_factor);
	jacobi_thread_expansion_pow = findLog_base2(thread_jacobi_expansion_factor);
	jacobi_stack_height = jacobi_block_expansion_pow + jacobi_thread_expansion_pow;
    jacobi_stack_height += 1; /*includes base with full size*/
	jacobi_scratch_stack_sizes = nullptr;
	jacobi_stack_numBlocks = nullptr;
	jacobi_stack_numThreads = nullptr;
	jacobi_stack_WH = nullptr;
    jacobi_scratch_stack = nullptr;
	jacobi_scratch = nullptr;
	b_stack = nullptr;
	Wx_stack = nullptr;
	Wy_stack = nullptr;
	jacobi_alpha = nullptr;
    jacobi_alpha_u = nullptr;
    jacobi_rbeta_u = nullptr;
	jacobi_delta_x = nullptr;
    jacobi_alpha_base = -delta_x * delta_x;
    jacobi_alpha_base_u = nu > 0.0 ? delta_x * delta_x / (nu * delta_t) : -1.0;
    jacobi_rbeta = 0.25;
    double jacobi_beta_u = 4 + jacobi_alpha_base_u;
    jacobi_rbeta_base_u = jacobi_beta_u > 0.0 ? 1.0 / jacobi_beta_u : -1.0;
    if (jacobi_stack_height > 0) {
		jacobi_scratch_stack_sizes = new int[jacobi_stack_height];
		jacobi_stack_numBlocks = new int[jacobi_stack_height];
		jacobi_stack_numThreads = new int[jacobi_stack_height];
		jacobi_stack_WH = new s_WH[jacobi_stack_height];
        jacobi_scratch_stack = new double* [jacobi_stack_height];
		b_stack = new double* [jacobi_stack_height];
        Wx_stack = new double* [jacobi_stack_height];
        Wy_stack = new double* [jacobi_stack_height];
        jacobi_alpha = new double[jacobi_stack_height];
        jacobi_alpha_u = new double[jacobi_stack_height];
        jacobi_rbeta_u = new double[jacobi_stack_height];
        jacobi_delta_x = new double[jacobi_stack_height];
        for (int i = 0; i < jacobi_stack_height; i++) {
            int reduction_factor = find_stack_WH_and_redFactor(jacobi_stack_WH[i], i, jacobi_stack_height);
            find_stack_BlocksNThreads(jacobi_stack_numBlocks[i], jacobi_stack_numThreads[i], reduction_factor);
			jacobi_delta_x[i] = delta_x * reduction_factor;
			jacobi_alpha[i] = -jacobi_delta_x[i] * jacobi_delta_x[i];
            jacobi_alpha_u[i] = nu > 0.0 ? jacobi_delta_x[i] * jacobi_delta_x[i] / (nu * delta_t) : -1.0;
            jacobi_beta_u = 4 + jacobi_alpha_u[i];
            jacobi_rbeta_u[i] = jacobi_beta_u > 0.0 ? 1.0 / jacobi_beta_u : -1.0;
        }
       
        for(int i=0; i<jacobi_stack_height; i++) {
			s_WH cur_wh = jacobi_stack_WH[i];
			jacobi_scratch_stack_sizes[i] = cur_wh.width * cur_wh.height;
			jacobi_scratch_stack[i] = nullptr;
			b_stack[i] = nullptr;
			Wx_stack[i] = nullptr;
			Wy_stack[i] = nullptr;
		}
    }
}
CUDAWrap::~CUDAWrap() {
    if(jacobi_delta_x != nullptr) {
        delete[] jacobi_delta_x;
	}
    if (jacobi_rbeta_u != nullptr) {
        delete[] jacobi_rbeta_u;
    }
    if (jacobi_alpha_u != nullptr) {
        delete[] jacobi_alpha_u;
    }
    if(jacobi_alpha != nullptr) {
        delete[] jacobi_alpha;
	}
    if (jacobi_scratch_stack != nullptr) {
		delete[] jacobi_scratch_stack;
    }
    if (jacobi_scratch_stack_sizes != nullptr) {
        delete[] jacobi_scratch_stack_sizes;
    }
    if (jacobi_stack_WH != nullptr) {
        delete[] jacobi_stack_WH;
    }
    if (jacobi_stack_numBlocks != nullptr) {
        delete[] jacobi_stack_numBlocks;
    }
    if (jacobi_stack_numThreads != nullptr) {
        delete[] jacobi_stack_numThreads;
    }
    if (b_stack != nullptr) {
        delete[] b_stack;
    }
    if (Wx_stack != nullptr) {
        delete[] Wx_stack;
    }
    if (Wy_stack != nullptr) {
        delete[] Wy_stack;
    }
}
int CUDAWrap::runCUDA(double* Ux, double* Uy, double* pressure, double* dye, int sim_frames, double jacobi_filter[], FluidAnimate* m_pfluid_animate) {

    unsigned int size = grid_width * grid_height;
    cudaError_t cudaStatus = cudaSuccess;
    if (initDevMem(size, jacobi_filter))
        cudaStatus = cudaSuccess;
    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaMalloc failed!");

    copyHostToDeviceMem(size, Ux, Uy, pressure, dye);

    int frames_run = 0;
    int frame_index = 0;
    int p_frame_index = 0;
    int p_advection_frame_index = 0;
	int dye_frame_index = 0;
    if (cudaStatus == cudaSuccess) {
        do {
			s_force force = m_pfluid_animate->updateForce(frames_run);
			s_force dye_brush = m_pfluid_animate->updateDye(frames_run);
            runFrame(frame_index, p_frame_index, p_advection_frame_index, dye_frame_index, force, dye_brush);
            frames_run++;
        } while (frames_run <= sim_frames);
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess)
            fprintf(stderr, "Cuda kernel launches failed:%s\n", cudaGetErrorString(cudaStatus));
    }
    add_full_grid(m_dev_p[p_advection_frame_index], m_dev_p0[p_frame_index]);
    copyDeviceToHostMem(size, Ux, Uy, pressure, dye, frame_index, p_advection_frame_index, dye_frame_index);

    releaseDevMem();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "runCUDA failed");
        return 1;
    }
    return 0;
}
void CUDAWrap::runFrame(
    int& frame_index,
    int& p_frame_index,
    int& p_advection_frame_index,
	int& dye_frame_index,
    s_force& force,
    s_force& dye_brush)
{
    int frame_in_index = frame_index;
    int p_frame_in_index = p_frame_index;
    int p_advection_frame_in_index = p_advection_frame_index;
    int dye_frame_in_index = dye_frame_index;
    advection(frame_in_index, dye_frame_in_index);
    reverseFrameIndex(frame_in_index);
    reverseFrameIndex(dye_frame_in_index);
    compute_pressure(m_dev_p, frame_in_index, p_advection_frame_in_index);
    reverseFrameIndex(p_advection_frame_in_index);
    subtract_pressure_gradient(m_dev_p, frame_in_index, p_advection_frame_in_index);
    reverseFrameIndex(frame_in_index);
    if (nu > 0.0) {
        viscous_diffusion(frame_in_index);
        reverseFrameIndex(frame_in_index);
    }
    if (force.active) {
        m_num_jacobi_loops = m_max_jacobi_force_loops;
        apply_force(frame_in_index, dye_frame_in_index, force, dye_brush);
        reverseFrameIndex(frame_in_index);
		reverseFrameIndex(dye_frame_in_index);
    }
    else {
        m_num_jacobi_loops = m_max_jacobi_loops;
    }
    compute_pressure(m_dev_p0, frame_in_index, p_frame_in_index);
    reverseFrameIndex(p_frame_in_index);
    subtract_pressure_gradient(m_dev_p0, frame_in_index, p_frame_in_index);
    reverseFrameIndex(frame_in_index);
    frame_index = frame_in_index;
    p_frame_index = p_frame_in_index;
    p_advection_frame_index = p_advection_frame_in_index;
	dye_frame_index = dye_frame_in_index;
}


void CUDAWrap::advection(int frame_index, int dye_frame_index) {
    cudaError_t cudaStatus = cudaSuccess;
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    s_frame_index frame_i = getFrameIndex(frame_index);
	s_frame_index dye_frame_i = getFrameIndex(dye_frame_index);
    advection_backtrace_to_indexes <<<numBlocks, numThreads>>>(
        advection_indexes[0],
        advection_indexes[1],
        advection_indexes[2],
        advection_indexes[3],
        advection_dist[0],
        advection_dist[1],
        advection_dist[2],
        advection_dist[3],
        m_dev_Ux[frame_i.in],
        m_dev_Uy[frame_i.in],
        delta_t,
        delta_x,
        grid_width,
        grid_height);
	cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess)
		fprintf(stderr, "advection_backtrace_to_indexes kernel launch failed:%s\n", cudaGetErrorString(cudaStatus));
	cudaStatus = cudaDeviceSynchronize();
	advection_Core<< <numBlocks_for_1D, numThreads_for_1D >> > (
        m_dev_Ux[frame_i.out], 
        m_dev_Uy[frame_i.out], 
		m_dev_dye[dye_frame_i.out],
        m_dev_Ux[frame_i.in], 
        m_dev_Uy[frame_i.in], 
		m_dev_dye[dye_frame_i.in],
        advection_indexes[0],
        advection_indexes[1],
        advection_indexes[2],
        advection_indexes[3],
        advection_dist[0],
        advection_dist[1],
        advection_dist[2],
        advection_dist[3]);
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "advection_Core kernel launch failed:%s\n", cudaGetErrorString(cudaStatus));
}
void CUDAWrap::viscous_diffusion(int frame_index) {
    s_frame_index frame_i = getFrameIndex(frame_index);
    jacobi_run(m_dev_Ux, m_dev_Ux[frame_index], nullptr, nullptr, frame_index);
    jacobi_run(m_dev_Uy, m_dev_Uy[frame_index], nullptr, nullptr, frame_index);
}
void CUDAWrap::viscous_diffusion_low(int frame_index) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    s_frame_index frame_i = getFrameIndex(frame_index);
	const double inv_delta_Xsqrd = 1.0 / (delta_x * delta_x);
    viscous_diffusion_core << <numBlocks, numThreads >> > (
        m_dev_Ux[frame_i.out],
        m_dev_Uy[frame_i.out],
        m_dev_Ux[frame_i.in],
        m_dev_Uy[frame_i.in],
        nu,
        inv_delta_Xsqrd,
        delta_t,
        grid_width,
		grid_height);
	cudaError_t cudaStatus = cudaDeviceSynchronize();
}
void CUDAWrap::apply_force(int frame_index, int dye_frame_index, s_force& force, s_force& dye_brush) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    s_frame_index frame_i = getFrameIndex(frame_index);
    s_frame_index dye_frame_i = getFrameIndex(dye_frame_index);
    double center_i = static_cast<double>(force.i);
    double center_j = static_cast<double>(force.j);
	double dye_center_i = static_cast<double>(dye_brush.i);
	double dye_center_j = static_cast<double>(dye_brush.j);
    double2 center = make_double2(center_i, center_j);
	double2 dye_center = make_double2(dye_center_i, dye_center_j);
    double2 Force = make_double2(force.Fx_c, force.Fy_c);
    applyForce_Core << <numBlocks, numThreads >> > (
        m_dev_Ux[frame_i.out],
        m_dev_Uy[frame_i.out],
		m_dev_dye[dye_frame_i.out],
        m_dev_Ux[frame_i.in],
        m_dev_Uy[frame_i.in],
		m_dev_dye[dye_frame_i.in],
        center,
        Force,
        dye_center,
        dye_brush.Fx_c,
        delta_t,
        force.inv_Rsqrd,
		dye_brush.inv_Rsqrd,
        grid_width);
    cudaError_t cudaStatus = cudaDeviceSynchronize();
}
void CUDAWrap::compute_pressure(double* p[], int frame_index, int p_frame_index) {
	divergence(m_dev_scratch, m_dev_Ux[frame_index], m_dev_Uy[frame_index]);
    /* jacobi with \alpha = -\deltax^2 and b = \frac{1}{\deltat} \nabla \cdot \vec{u} and \beta = 4 */
    //alpha = -delta_x * delta_x;
    //rbeta = 0.25;
    jacobi_run(p, m_dev_scratch, m_dev_Ux[frame_index], m_dev_Uy[frame_index], p_frame_index);
}
void CUDAWrap::subtract_pressure_gradient(double* p[], int frame_index, int p_frame_index) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    s_frame_index frame_i = getFrameIndex(frame_index);
    s_frame_index p_frame_i = getFrameIndex(p_frame_index);
    static double inv_2delta_x = 1.0 / (2.0 * delta_x);
    subtractGradient << <numBlocks, numThreads >> > (
        m_dev_Ux[frame_i.out],
        m_dev_Uy[frame_i.out],
        m_dev_Ux[frame_i.in],
        m_dev_Uy[frame_i.in],
        p[p_frame_i.in],
        inv_2delta_x,
        grid_width,
        grid_height);
    cudaError_t cudaStatus = cudaDeviceSynchronize();
}

void CUDAWrap::jacobi_fill_stacks(const double* frame_in, const double* b, const double* Wx, const double* Wy) {
	cudaError_t cudaStatus = cudaSuccess;
    int jacobi_stack_max_index = jacobi_stack_height - 1;
    copy_memory << <numBlocks_for_1D, numThreads_for_1D >> > (jacobi_scratch_stack[jacobi_stack_max_index], frame_in);
    copy_memory << <numBlocks_for_1D, numThreads_for_1D >> > (b_stack[jacobi_stack_max_index], b);
    if (Wx != nullptr && Wy != nullptr) {
        copy_memory << <numBlocks_for_1D, numThreads_for_1D >> > (Wx_stack[jacobi_stack_max_index], Wx);
        copy_memory << <numBlocks_for_1D, numThreads_for_1D >> > (Wy_stack[jacobi_stack_max_index], Wy);
    }
    int r_grid_width = grid_width;
    int r_grid_height = grid_height;
    int base_grid_width = grid_width;
    int base_grid_height = grid_height;
    for (int i = 0; i < jacobi_stack_max_index; i++)
    {
        r_grid_width /= 2;
        r_grid_height /= 2;
        int stack_base_i = jacobi_stack_max_index - i;
        int r_stack_i = jacobi_stack_max_index - i - 1;
        int numBlocks_s = jacobi_stack_numBlocks[r_stack_i];
        int numThreads_s = jacobi_stack_numThreads[r_stack_i];
        dim3 numBlocks(numBlocks_s, numBlocks_s);
        dim3 numThreads(numThreads_s, numThreads_s);
        if (Wx != nullptr && Wy != nullptr) {
            Xgrid_reduction_x4 << <numBlocks, numThreads >> > (
                jacobi_scratch_stack[r_stack_i],
                b_stack[r_stack_i],
                Wx_stack[r_stack_i],
                Wy_stack[r_stack_i],
                r_grid_width,
                r_grid_height,
                jacobi_scratch_stack[stack_base_i],
                b_stack[stack_base_i],
                Wx_stack[stack_base_i],
                Wy_stack[stack_base_i],
                base_grid_width,
                base_grid_height,
                m_filter);
        }
        else {
            Xgrid_reduction_dup << <numBlocks, numThreads >> > (
                jacobi_scratch_stack[r_stack_i],
                b_stack[r_stack_i],
                r_grid_width,
                r_grid_height,
                jacobi_scratch_stack[stack_base_i],
                base_grid_width,
                base_grid_height,
                m_filter);
        }
		cudaStatus = cudaGetLastError();
        if(cudaStatus!= cudaSuccess)
			fprintf(stderr, "Xgrid_reduction_x4 kernel launch failed:%s\n", cudaGetErrorString(cudaStatus));
        cudaStatus = cudaDeviceSynchronize();
        base_grid_width /= 2;
        base_grid_height /= 2;
    }
}

void CUDAWrap::jacobi_send_frame_down_stack(const double* hi_frame, int hi_frame_width_height, int lo_stack_index) {
	cudaError_t cudaStatus = cudaSuccess;
    int numBlocks_s = jacobi_stack_numBlocks[lo_stack_index];
    int numThreads_s = jacobi_stack_numThreads[lo_stack_index];
    dim3 numBlocks(numBlocks_s, numBlocks_s);
    dim3 numThreads(numThreads_s, numThreads_s);
    int lo_grid_width_height = numThreads_s * numBlocks_s;
    Xgrid_expansion << <numBlocks, numThreads >> > (
        jacobi_scratch_stack[lo_stack_index],
        lo_grid_width_height,
        lo_grid_width_height,
        hi_frame,
        hi_frame_width_height,
        hi_frame_width_height);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "Xgrid_expansion kernel launch failed:%s\n", cudaGetErrorString(cudaStatus));
    cudaStatus = cudaDeviceSynchronize();
}

void CUDAWrap::jacobi_frame(
    double* frame_out,
    const double* frame_in,
    const double* b,
    const double* Wx,
    const double* Wy,
    const double& alpha,
    const double& rbeta,
    const int& numBlocks_s,
    const int& numThreads_s)
{
    cudaError_t cudaStatus = cudaSuccess;
    dim3 numBlocks(numBlocks_s, numBlocks_s);
    dim3 numThreads(numThreads_s, numThreads_s);
    int frame_wh = numBlocks_s * numThreads_s;
    jacobi << <numBlocks, numThreads >> > (frame_out, frame_in, b, alpha, rbeta, frame_wh, frame_wh);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "jacobi kernel launch failed:%s\n", cudaGetErrorString(cudaStatus));
    cudaStatus = cudaDeviceSynchronize();
    if (Wx != nullptr && Wy != nullptr) {
        jacobi_boundary_pressure << <numBlocks_side, numThreads_side >> > (frame_out, Wx, Wy, frame_in, b, alpha, rbeta, delta_x, frame_wh, frame_wh);
        cudaStatus = cudaDeviceSynchronize();
    }
}

void CUDAWrap::jacobi_loop(
    double* X[],
    const double* b,
    double alpha,
    double rbeta,
    int numBlocks_s,
    int numThreads_s,
    const double* Wx,
    const double* Wy)
{
    s_frame_index frame_i = getFrameIndex(0);
    int jacobi_loop_index = 0;
    do {
        jacobi_frame(X[frame_i.out], X[frame_i.in], b, Wx, Wy, alpha, rbeta, numBlocks_s, numThreads_s);
        swapFrameIndexes(frame_i);/*results are now at frame_i.in */
        jacobi_loop_index++;
    } while (jacobi_loop_index <= m_num_jacobi_loops);
    fixFramePointers(X, frame_i, 0);/* results are now in X[0] */
}
void CUDAWrap::jacobi_run(
    double* X[],
    const double* b,
    const double* Wx,
    const double* Wy,
    int frame_index)
{
    s_frame_index frame_i = getFrameIndex(frame_index);
    int original_frame_in_index = frame_i.in;
    if (jacobi_stack_height > 0) {
        jacobi_fill_stacks(X[frame_i.in], b, Wx, Wy);
        int jacobi_stack_max_index = jacobi_stack_height - 1;
        for (int i_stack = 0; i_stack < jacobi_stack_max_index; i_stack++) {
            double* frame_swap_ptrs[] = { jacobi_scratch_stack[i_stack], jacobi_scratch };
            int numBlocks_s = jacobi_stack_numBlocks[i_stack];
            int numThreads_s = jacobi_stack_numThreads[i_stack];
            double alpha, rbeta;
            if (Wx != nullptr && Wy != nullptr) {
                alpha = jacobi_alpha[i_stack];
                rbeta = jacobi_rbeta;
                jacobi_loop(frame_swap_ptrs, b_stack[i_stack], alpha, rbeta, numBlocks_s, numThreads_s, Wx_stack[i_stack], Wy_stack[i_stack]);
            }
            else {
                alpha = jacobi_alpha_u[i_stack];
                rbeta = jacobi_rbeta_u[i_stack];
                jacobi_loop(frame_swap_ptrs, b_stack[i_stack], alpha, rbeta, numBlocks_s, numThreads_s, nullptr, nullptr);
            }
            s_WH wh = jacobi_stack_WH[i_stack];
            int frame_width_height = wh.width;
            jacobi_send_frame_down_stack(frame_swap_ptrs[0], frame_width_height, i_stack + 1);
        }
        copyMemory_for_standard_grid(X[0], jacobi_scratch_stack[jacobi_stack_max_index]);
    }
    double _alpha_base, _rbeta_base; 
    if (Wx != nullptr && Wy != nullptr) {
        _alpha_base = jacobi_alpha_base;
        _rbeta_base = jacobi_rbeta;
    }
    else
    {
        _alpha_base = jacobi_alpha_base_u;
        _rbeta_base = jacobi_rbeta_base_u;
    }
    jacobi_loop(X, b, _alpha_base, _rbeta_base, numBlocks_side, numThreads_side, Wx, Wy);/* results are in X[0]*/
    frame_i.in = 1;
    frame_i.out = 0;
    fixFramePointers(X, frame_i, original_frame_in_index);/*results are in original frame out index */
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
void CUDAWrap::bilinearAprox_scaledFrame(double* Ux_scaled, double* Uy_scaled, double* dye_scaled, const double* Ux, const double* Uy, const double* dye, int scale_factor) {
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

void CUDAWrap::gradient(double* dp_dx, double* dp_dy, double* p[], int p_frame_index) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
	s_frame_index p_frame_i = getFrameIndex(p_frame_index);
    static double inv_2delta_x = 1.0 / (2.0*delta_x);
    gradient_core << <numBlocks, numThreads >> > (
        dp_dx,
        dp_dy,
        p[p_frame_i.in],
        inv_2delta_x,
        grid_width,
        grid_height);
	cudaError_t cudaStatus = cudaDeviceSynchronize();
}
void CUDAWrap::divergence(double* div, const double* Ux, const double* Uy) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    static double inv_2delta_x = 1.0 / (2.0 * delta_x);
    divergence_Core << <numBlocks, numThreads >> > (div, Ux, Uy, inv_2delta_x, grid_width, grid_height);
    cudaError_t cudaStatus = cudaDeviceSynchronize();
}
void CUDAWrap::derivative(double* dX, double* dY, const double* Ux, const double* Uy) {
    dim3 numBlocks(numBlocks_side, numBlocks_side);
    dim3 numThreads(numThreads_side, numThreads_side);
    static double inv_2delta_x = 1.0 / (2.0 * delta_x);
    derivative_Core << <numBlocks, numThreads >> > (dX, dY, Ux, Uy, inv_2delta_x, grid_width, grid_height);
	cudaError_t cudaStatus = cudaDeviceSynchronize();
}
void CUDAWrap::laplacian(double* lap, double* p[], double* scratch_x, double* scratch_y, int p_frame_index) {
    gradient(scratch_x, scratch_y, p, p_frame_index);
	dim3 numBlocks(numBlocks_side, numBlocks_side);
	dim3 numThreads(numThreads_side, numThreads_side);
    static double inv_2delta_x = 1.0 / (2.0 * delta_x);
	divergence_Core << <numBlocks, numThreads >> > (lap, scratch_x, scratch_y, inv_2delta_x, grid_width, grid_height);
	cudaError_t cudaStatus = cudaDeviceSynchronize();
}



bool CUDAWrap::initDevMem(int size, double jacobi_filter[]) {
    cudaError_t cudaStatus = cudaSuccess;
    for (int i = 0; i < 2; i++) {
        m_dev_Ux[i] = 0;
        m_dev_Uy[i] = 0;
        m_dev_p[i] = 0;
        m_dev_p0[i] = 0;
    }
    m_dev_scratch = 0;
    for(int i=0; i<4; i++){
        advection_indexes[i] = 0;
		advection_dist[i] = 0;
	}

    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaSetDevice failed!");

    for (int coord_index = 0; coord_index < 2; coord_index++) {
        if (cudaStatus == cudaSuccess)
            cudaStatus = cudaMalloc((void**)&m_dev_Ux[coord_index], size * sizeof(double));
        if (cudaStatus == cudaSuccess)
            cudaStatus = cudaMalloc((void**)&m_dev_Uy[coord_index], size * sizeof(double));
        if (cudaStatus == cudaSuccess)
			cudaStatus = cudaMalloc((void**)&m_dev_p0[coord_index], size * sizeof(double));
        if (cudaStatus == cudaSuccess)
            cudaStatus = cudaMalloc((void**)&m_dev_p[coord_index], size * sizeof(double));
        if (cudaStatus == cudaSuccess)
            cudaStatus = cudaMalloc((void**)&m_dev_dye[coord_index], size * sizeof(double));
    }
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&m_dev_scratch, size * sizeof(double));
    for(int coord_index=0; coord_index < 4; coord_index++) {
        if (cudaStatus == cudaSuccess)
            cudaStatus = cudaMalloc((void**)&advection_indexes[coord_index], size * sizeof(int));
        if (cudaStatus == cudaSuccess)
            cudaStatus = cudaMalloc((void**)&advection_dist[coord_index], size * sizeof(double));
	}
    for (int i = 0; i < jacobi_stack_height; i++) {
        jacobi_scratch_stack[i] = 0;
        b_stack[i] = 0;
        Wx_stack[i] = 0;
        Wy_stack[i] = 0;
        if (cudaStatus == cudaSuccess)
            cudaStatus = cudaMalloc((void**)&jacobi_scratch_stack[i], jacobi_scratch_stack_sizes[i] * sizeof(double));
        if (cudaStatus == cudaSuccess)
            cudaStatus = cudaMalloc((void**)&b_stack[i], jacobi_scratch_stack_sizes[i] * sizeof(double));
        if (cudaStatus == cudaSuccess)
            cudaStatus = cudaMalloc((void**)&Wx_stack[i], jacobi_scratch_stack_sizes[i] * sizeof(double));
        if (cudaStatus == cudaSuccess)
            cudaStatus = cudaMalloc((void**)&Wy_stack[i], jacobi_scratch_stack_sizes[i] * sizeof(double));
    }
    jacobi_scratch = 0;
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&jacobi_scratch, size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&m_filter, g_jacobi_filter_size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(m_filter, jacobi_filter, g_jacobi_filter_size * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess)
        return false;
    return true;
}
bool CUDAWrap::copyHostToDeviceMem(int size, double* Ux, double* Uy, double* pressure, double* dye) {
    cudaError_t cudaStatus = cudaSuccess;
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(m_dev_Ux[0], Ux, size * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(m_dev_Uy[0], Uy, size * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(m_dev_p0[0], pressure, size * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(m_dev_p[0], pressure, size * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus == cudaSuccess)
		cudaStatus = cudaMemcpy(m_dev_dye[0], dye, size * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaMemcpy of host to device failed!");
    return cudaStatus == cudaSuccess;
}
bool CUDAWrap::copyDeviceToHostMem(int size, double* Ux, double* Uy, double* pressure, double* dye, int frame_index, int p_frame_index, int dye_frame_index) {
    cudaError_t cudaStatus = cudaSuccess;

    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(Ux, m_dev_Ux[frame_index], size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(Uy, m_dev_Uy[frame_index], size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(pressure, m_dev_p[p_frame_index], size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(dye, m_dev_dye[dye_frame_index], size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaMemcpy of device to host failed!");
    return cudaStatus == cudaSuccess;
}
void CUDAWrap::releaseDevMem() {
    for (int i = 0; i < 4; i++) {
        cudaFree(advection_indexes[i]);
        cudaFree(advection_dist[i]);
    }
    cudaFree(jacobi_scratch);
    for (int i = 0; i < jacobi_stack_height; i++) {
        cudaFree(jacobi_scratch_stack[i]);
        cudaFree(b_stack[i]);
        cudaFree(Wx_stack[i]);
        cudaFree(Wy_stack[i]);
    }
    for (int coord_index = 0; coord_index < 2; coord_index++) {
        cudaFree(m_dev_Ux[coord_index]);
        cudaFree(m_dev_Uy[coord_index]);
        cudaFree(m_dev_p[coord_index]);
		cudaFree(m_dev_p0[coord_index]);
		cudaFree(m_dev_dye[coord_index]);
    }
    cudaFree(m_dev_scratch);
}

void CUDAWrap::add_full_grid(double* dest, const double* src) {
    add_values << <numBlocks_for_1D, numThreads_for_1D >> > (dest, src);
    cudaError_t cudaStatus = cudaDeviceSynchronize();
}
void CUDAWrap::copyMemory_for_standard_grid(double* copy, const double* orig) {
    copy_memory << <numBlocks_for_1D, numThreads_for_1D >> > (copy, orig);
    cudaError_t cudaStatus = cudaDeviceSynchronize();
}
int CUDAWrap::findLog_base2(int val) {
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
int CUDAWrap::find2Pow(int pow) {
    int val = 1;
    for (int i = 0; i < pow; i++)
        val *= 2;
    return val;
}
int CUDAWrap::find_stack_WH_and_redFactor(s_WH& wh, int stack_index, int stack_height) {
	int pow2 = (stack_height-1) - stack_index;
	int reduction_factor = find2Pow(pow2);
    wh.width = grid_width / reduction_factor;
    wh.height = grid_height / reduction_factor;
    return reduction_factor;
}
void CUDAWrap::find_stack_BlocksNThreads(int& numBlocks_s, int& numThreads_s, int reduction_factor) {
    numBlocks_s = numBlocks_side;
	numThreads_s = numThreads_side;
    int red_fac_done = 1;
    while (numBlocks_s > jacobi_minBlocks_side && red_fac_done < reduction_factor) {
        numBlocks_s /= 2;
		red_fac_done *= 2;
    }
    while(red_fac_done<reduction_factor) {
        numThreads_s /= 2;
        red_fac_done *= 2;
	}
}
