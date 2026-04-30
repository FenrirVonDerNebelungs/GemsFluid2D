#include "Test.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

Test::Test(int num_display_frames, int blow_factor, int loop_modulus_divider) :
    m_current_frame(0),
    m_Ux_max(0.0),
    m_Uy_max(0.0),
    m_U_max(0.0),
    m_p_max(0.0),
    m_mouse_clicks(0),
    m_blow_factor(blow_factor),
    m_loop_modulus_divider(loop_modulus_divider),
    m_pCurrent_GenImage(nullptr),
    m_current_max(0.0)
{
    m_pCUDA_wrap = new CUDAWrap();
	int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
	m_Ux = new double[size];
	m_Uy = new double[size];
	m_p = new double[size];
    int blown_size = size * m_blow_factor * m_blow_factor;
    m_Ux_bilinear = new double[blown_size];
    m_Uy_bilinear = new double[blown_size];
    m_pPyTrans = new PyTrans();
	m_pFluid_animate = new FluidAnimate();
	s_WH wh = getGridWidthHeight();
	m_pFluid_animate->init(wh);
	m_force = m_pFluid_animate->getForce(0);
    for (int i = 0; i < 4; i++) {
		m_host_scratch[i] = new double[size];
    }
    m_pDraw = new drawTest(m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height);
}
Test::~Test() {
    if (m_pDraw != nullptr)
        delete m_pDraw;
    for (int i = 0; i < 4; i++) {
        if (m_host_scratch[i] != nullptr)
            delete[] m_host_scratch[i];
    }
    if (m_pFluid_animate != nullptr)
        delete m_pFluid_animate;
    if (m_pPyTrans != nullptr)
        delete m_pPyTrans;
    if (m_Uy_bilinear != nullptr)
        delete [] m_Uy_bilinear;
    if (m_Ux_bilinear != nullptr)
        delete[] m_Ux_bilinear;

    if(m_p!= nullptr)
		delete[] m_p;
    if(m_Uy!=nullptr)
		delete[] m_Uy;
	if (m_Ux != nullptr)
		delete[] m_Ux;
	if (m_pCUDA_wrap != nullptr)
        delete m_pCUDA_wrap;
}
int Test::runTest(int sim_frames) {
    int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    std::memset(m_Ux, 0, size * sizeof(double));
    std::memset(m_Uy, 0, size * sizeof(double));
    std::memset(m_p, 0, size * sizeof(double));
    runCUDA(m_Ux, m_Uy, m_p, m_force, sim_frames);
    return 0;
}
int Test::runCUDA(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames) {
    unsigned int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    cudaError_t cudaStatus = cudaSuccess;
    if (m_pCUDA_wrap->initDevMem(size))
        cudaStatus = cudaSuccess;
    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaMalloc failed!");


    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(m_pCUDA_wrap->m_dev_Ux[0], Ux, size * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(m_pCUDA_wrap->m_dev_Uy[0], Uy, size * sizeof(double), cudaMemcpyHostToDevice);
    int frames_run = 0;
    int frame_index = 0;
    int p_frame_index = 0;
    if (cudaStatus == cudaSuccess) {
        int cache_full_len = getCacheLen();
        m_pPyTrans->init("../Dat/frames.dat", cache_full_len, 0);
        do {
            runFrame(
                m_pCUDA_wrap->m_dev_Ux, 
                m_pCUDA_wrap->m_dev_Uy, 
                m_pCUDA_wrap->m_dev_p, 
                m_pCUDA_wrap->m_dev_scratch, 
                frame_index, 
                p_frame_index, 
                force);
            m_pPyTrans->resetCache();//AndWrite();
            frames_run++;
        } while (frames_run <= sim_frames);
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess)
            fprintf(stderr, "Cuda kernel launches failed:%s\n", cudaGetErrorString(cudaStatus));
        m_pPyTrans->release();
    }
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(Ux, m_pCUDA_wrap->m_dev_Ux[frame_index], size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(Uy, m_pCUDA_wrap->m_dev_Uy[frame_index], size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(pressure, m_pCUDA_wrap->m_dev_p[p_frame_index], size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaMemcpy failed!");

    m_pCUDA_wrap->releaseDevMem();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "runCUDA failed");
        return 1;
    }
    return 0;
}
void Test::runFrame(double* Ux[], double* Uy[], double* p[], double* scratch, int& frame_index, int& p_frame_index, s_force& force) {
    int frame_in_index = frame_index;
    int p_frame_in_index = p_frame_index;
	send_frame_to_host(m_Ux, Ux[frame_index]);
	send_frame_to_host(m_Uy, Uy[frame_index]);
	send_frame_to_host(m_p, p[p_frame_index]);
    m_pPyTrans->cacheGrid(m_Ux, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::U_code, n_PyTrans::X_code, n_PyTrans::start_frame_code);
    m_pPyTrans->cacheGrid(m_Uy, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::U_code, n_PyTrans::Y_code, n_PyTrans::start_frame_code);
    m_pCUDA_wrap->advection(Ux, Uy, frame_in_index);
    reverseFrameIndex(frame_in_index);
    send_frame_to_host(m_Ux, Ux[frame_in_index]);
    send_frame_to_host(m_Uy, Uy[frame_in_index]);
    m_pPyTrans->cacheGrid(m_Ux, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::W_code, n_PyTrans::X_code, n_PyTrans::after_advection_code);
    m_pPyTrans->cacheGrid(m_Uy, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::W_code, n_PyTrans::Y_code, n_PyTrans::after_advection_code);
    m_pPyTrans->cacheGrid(m_p, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::P_code, n_PyTrans::Scalar_code, n_PyTrans::start_frame_code);
    m_pCUDA_wrap->compute_pressure(Ux, Uy, p, scratch, frame_in_index, p_frame_in_index);
    reverseFrameIndex(p_frame_in_index);
	send_frame_to_host(m_p, p[p_frame_in_index]);
    m_pPyTrans->cacheGrid(m_p, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::P_code, n_PyTrans::Scalar_code, n_PyTrans::after_advection_code);
    m_pCUDA_wrap->subtract_pressure_gradient(Ux, Uy, p, frame_in_index, p_frame_in_index);
    reverseFrameIndex(frame_in_index);
    send_frame_to_host(m_Ux, Ux[frame_index]);
    send_frame_to_host(m_Uy, Uy[frame_index]);
    m_pPyTrans->cacheGrid(m_Ux, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::U_code, n_PyTrans::X_code, n_PyTrans::after_advection_code);
    m_pPyTrans->cacheGrid(m_Uy, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::U_code, n_PyTrans::Y_code, n_PyTrans::after_advection_code);
    //m_pCUDA_wrap->viscous_diffusion(Ux, Uy, frame_in_index);
    //reverseFrameIndex(frame_in_index);
    if (force.active) {
        m_pCUDA_wrap->apply_force(Ux, Uy, frame_in_index, force);
        reverseFrameIndex(frame_in_index);
    }
    send_frame_to_host(m_Ux, Ux[frame_index]);
    send_frame_to_host(m_Uy, Uy[frame_index]);
    m_pPyTrans->cacheGrid(m_Ux, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::W_code, n_PyTrans::X_code, n_PyTrans::after_force_code);
    m_pPyTrans->cacheGrid(m_Uy, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::W_code, n_PyTrans::Y_code, n_PyTrans::after_force_code);
    compute_pressure(Ux, Uy, p, scratch, frame_in_index, p_frame_in_index);
    reverseFrameIndex(p_frame_in_index);
	s_frame_index test_frame_i = getFrameIndex(frame_in_index);
	m_pCUDA_wrap->gradient(Ux[test_frame_i.out], Uy[test_frame_i.out], p, frame_in_index);
	send_frame_to_host(m_Ux, Ux[frame_index]);
	send_frame_to_host(m_Uy, Uy[frame_index]);
	m_pPyTrans->cacheGrid(m_Ux, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::gradP_code, n_PyTrans::X_code, n_PyTrans::after_pressure_code);
	m_pPyTrans->cacheGrid(m_Uy, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::gradP_code, n_PyTrans::Y_code, n_PyTrans::after_pressure_code);
    m_pCUDA_wrap->subtract_pressure_gradient(Ux, Uy, p, frame_in_index, p_frame_in_index);
    reverseFrameIndex(frame_in_index);
    send_frame_to_host(m_Ux, Ux[frame_index]);
    send_frame_to_host(m_Uy, Uy[frame_index]);
    m_pPyTrans->cacheGrid(m_Ux, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::U_code, n_PyTrans::X_code, n_PyTrans::end_frame_code);
    m_pPyTrans->cacheGrid(m_Uy, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::U_code, n_PyTrans::Y_code, n_PyTrans::end_frame_code);
    /*14 cachegrid ops of (grid_width*grid_height + header each) */
    frame_index = frame_in_index;
    p_frame_index = p_frame_in_index;
}
void Test::compute_pressure(double* Ux[], double* Uy[], double* p[], double* scratch, int frame_index, int p_frame_index) {
    m_pCUDA_wrap->divergence(scratch, Ux[frame_index], Uy[frame_index]);
	send_frame_to_host(m_p, scratch);
	m_pPyTrans->cacheGrid(m_p, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::DivW_code, n_PyTrans::Scalar_code, n_PyTrans::after_force_code);
    /* jacobi with \alpha = -\deltax^2 and b = \frac{1}{\deltat} \nabla \cdot \vec{u} and \beta = 4 */
    static double alpha = -m_pCUDA_wrap->delta_x * m_pCUDA_wrap->delta_x;
    static double rbeta = 0.25;
    jacobi_run(p, scratch, Ux[frame_index], Uy[frame_index], p_frame_index, alpha, rbeta);
    /* 1 cache grid ops grid_width*grid_height + header  */
}
void Test::jacobi_loop(
    double* X[],
    const double* b,
    double alpha,
    double rbeta,
    int numBlocks_s,
    int numThreads_s,
    const double* Wx,
    const double* Wy,
    int stack_index)
{
    s_frame_index frame_i = getFrameIndex(0);
    int num_jacobi_loops = 1;
    do {
        m_pCUDA_wrap->jacobi_frame(X[frame_i.out], X[frame_i.in], b, Wx, Wy, alpha, rbeta, numBlocks_s, numThreads_s);
        swapFrameIndexes(frame_i);/*results are now at frame_i.in */
        if ((num_jacobi_loops-1) % m_loop_modulus_divider == 0) {
            int frame_width = numBlocks_s * numThreads_s;
            int frame_height = frame_width;
            int frame_len = frame_width * frame_height;
            send_frame_to_host(m_p, X[frame_i.in], frame_len);
            m_pPyTrans->cacheGrid(
                m_p,
                frame_width,
                frame_height,
                m_current_frame,
                n_PyTrans::P_code,
                n_PyTrans::Scalar_code,
                n_PyTrans::jacobi_loop_frame_code,
                num_jacobi_loops,
                stack_index);
        }
        num_jacobi_loops++;
    } while (num_jacobi_loops < m_pCUDA_wrap->max_jacobi_loops);
    fixFramePointers(X, frame_i, 0);/* results are now in X[0] */
    /** (max_jacobi_loops/loop_modulus_divider)*(header + cur_frame_len) **/
}
void Test::jacobi_run(
    double* X[],
    const double* b,
    const double* Wx,
    const double* Wy,
    int frame_index,
    const double& alpha,
    const double& rbeta)
{
    s_frame_index frame_i = getFrameIndex(frame_index);
    int original_frame_in_index = frame_i.in;
    if (m_pCUDA_wrap->jacobi_stack_height > 0) {
        m_pCUDA_wrap->jacobi_fill_stacks(X[frame_i.in], b, Wx, Wy);
		send_jacobi_stacks_to_cache(0);
        int jacobi_stack_max_index = m_pCUDA_wrap->jacobi_stack_height - 1;
        for (int i_stack = 0; i_stack < jacobi_stack_max_index; i_stack++) {
            double* frame_swap_ptrs[] = { m_pCUDA_wrap->jacobi_scratch_stack[i_stack], m_pCUDA_wrap->jacobi_scratch };
            int numBlocks_s = m_pCUDA_wrap->jacobi_stack_numBlocks[i_stack];
            int numThreads_s = m_pCUDA_wrap->jacobi_stack_numThreads[i_stack];
            jacobi_loop(frame_swap_ptrs, m_pCUDA_wrap->b_stack[i_stack], alpha, rbeta, numBlocks_s, numThreads_s, m_pCUDA_wrap->Wx_stack[i_stack], m_pCUDA_wrap->Wy_stack[i_stack], i_stack);
			s_WH wh = m_pCUDA_wrap->jacobi_stack_WH[i_stack];
            send_frame_to_host(m_p, frame_swap_ptrs[0], wh.width * wh.height);
            m_pPyTrans->cacheGrid(
                m_p,
                wh.width,
                wh.height,
                m_current_frame,
                n_PyTrans::P_code,
                n_PyTrans::Scalar_code,
                n_PyTrans::jacobi_frame_code,
                m_pCUDA_wrap->max_jacobi_loops - 1,
                i_stack);
            int frame_width_height = numBlocks_s * numThreads_s;
            m_pCUDA_wrap->jacobi_send_frame_down_stack(frame_swap_ptrs[0], frame_width_height, i_stack + 1);
            /* (jacobi_stack_height-1) * (varying stack width*heights)*header_len **/
        }
        m_pCUDA_wrap->copyMemory_for_standard_grid(X[0], m_pCUDA_wrap->jacobi_scratch_stack[jacobi_stack_max_index]);
    }
    jacobi_loop(X, b, alpha, rbeta, m_pCUDA_wrap->numBlocks_side, m_pCUDA_wrap->numThreads_side, Wx, Wy,0);/* results are in X[0]*/
    /* jacobi_stack_height*(max_jacobi_loops/loop_modulus_divider)*(header + full_stack_len) */
	send_frame_to_host(m_p, X[0]);
    m_pPyTrans->cacheGrid(
        m_p,
        m_pCUDA_wrap->grid_width,
        m_pCUDA_wrap->grid_height,
        m_current_frame,
        n_PyTrans::P_code,
        n_PyTrans::Scalar_code,
        n_PyTrans::jacobi_frame_code,
        m_pCUDA_wrap->max_jacobi_loops - 1,
        0);
    /* jacobi_stack_height*header_len + full jacobi stack len */
    frame_i.in = 1;
    frame_i.out = 0;
    fixFramePointers(X, frame_i, original_frame_in_index);/*results are in original frame out index */
 
    /** total for this function and sub functions **/
    /* (5 + max_jacobi_loops/loop_modulus_divider)*jacobi_stack_height*(header len + full jacobi stack len) */
}
void Test::send_jacobi_stacks_to_cache(int jacobi_frame) {
    for (int i_stack = 0; i_stack < m_pCUDA_wrap->jacobi_stack_height; i_stack++) {
		s_WH wh = m_pCUDA_wrap->jacobi_stack_WH[i_stack];
        int frame_width_height = wh.width * wh.height;
        send_frame_to_host(m_host_scratch[0], m_pCUDA_wrap->jacobi_scratch_stack[i_stack], frame_width_height);
		send_frame_to_host(m_host_scratch[1], m_pCUDA_wrap->b_stack[i_stack], frame_width_height);
		send_frame_to_host(m_host_scratch[2], m_pCUDA_wrap->Wx_stack[i_stack], frame_width_height);
		send_frame_to_host(m_host_scratch[3], m_pCUDA_wrap->Wy_stack[i_stack], frame_width_height);
        m_pPyTrans->cacheGrid(
            m_host_scratch[0],
            wh.width,
            wh.height,
            m_current_frame,
            n_PyTrans::P_code,
            n_PyTrans::Scalar_code,
            n_PyTrans::jacobi_stack_fill_code,
            jacobi_frame,
            i_stack
        );
        m_pPyTrans->cacheGrid(
            m_host_scratch[1],
            wh.width,
            wh.height,
            m_current_frame,
            n_PyTrans::DivW_code,
            n_PyTrans::Scalar_code,
            n_PyTrans::jacobi_stack_fill_code,
            jacobi_frame,
            i_stack
        );
        m_pPyTrans->cacheGrid(
            m_host_scratch[2],
            wh.width,
            wh.height,
            m_current_frame,
            n_PyTrans::W_code,
            n_PyTrans::X_code,
            n_PyTrans::jacobi_stack_fill_code,
            jacobi_frame,
            i_stack
        );
        m_pPyTrans->cacheGrid(
            m_host_scratch[3],
            wh.width,
            wh.height,
            m_current_frame,
            n_PyTrans::W_code,
            n_PyTrans::Y_code,
            n_PyTrans::jacobi_stack_fill_code,
            jacobi_frame,
            i_stack
        );
	}
    /* 4 * (header len + full jacobi_stack len) */
}
int Test::getCacheLen() {
    /* 14 full size from runFrame */
       /*compute pressure-> 1 full size frame*/
          /*jacobi run -> jacobi_stack_height * header + full_jacobi_stack_len */
               /* send_jacobi_stacks_to_cache -> 4*(jacobi_stack_height * header + full_jacobi_stack_len)*/
               /*jacobi_loop (max_jacobi_loops/loop_modulus_divider)*(header + stack_len0 + header + stack_len1 + ... )
                             (max_jacobi_loops/loop_modulus_divider)*jacobi_stack_height*(header) + (max_jacobi_loops/loop_modulus_divider)*full_jacobi_stack_len */
    /* full_jacobi_stack len  = min_stack_len + 4*min_stack_len + 4^2min_stack_len + 4^3min_stack_len ... largest_stack_len  
       largest_stack_len = full_size_frame 
       4^N = largest_stack_len/min_stack_len 
       full_jacobi_stack_len = min_stack_len*(1 + 4^2 +4^3 +.. 4^N)
       S_N = 1+r^2+r^3+..r^N = (1-r^{N+1})/(1-r)
       full_jacobi_stack_len = min_stack_len* \frac{1-(4^{N})*4}{-3} */
	int full_size_frame_len = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height + n_PyTrans::obj_header_len;
	int cache_len = 14 * full_size_frame_len; /* dumps directly from runframe */
    cache_len += full_size_frame_len; /*1 frame direclty in compute pressure*/

    int min_sized_stack_frame_len_data_only = m_pCUDA_wrap->jacobi_scratch_stack_sizes[0];
	int _4_pow_N = (m_pCUDA_wrap->jacobi_scratch_stack_sizes[m_pCUDA_wrap->jacobi_stack_height - 1]) / min_sized_stack_frame_len_data_only;
	int full_stack_frame_sum_len_data_only = min_sized_stack_frame_len_data_only * (1 - _4_pow_N * 4) / -3;

    cache_len += 4 * (m_pCUDA_wrap->jacobi_stack_height * n_PyTrans::obj_header_len + full_stack_frame_sum_len_data_only);
    int num_cache_writes_per_jacobi_loop = m_pCUDA_wrap->max_jacobi_loops / m_loop_modulus_divider;
    cache_len += num_cache_writes_per_jacobi_loop * (m_pCUDA_wrap->jacobi_stack_height * n_PyTrans::obj_header_len + full_stack_frame_sum_len_data_only);
    
    return cache_len;
}
double* Test::getFrameToDisplay() {
    m_pCurrent_GenImage = m_pDraw->getGenImage();
    int mouse_clicks = m_mouse_clicks;
    m_mouse_clicks++;
    switch (mouse_clicks) {
    case 0:
        m_pDraw->setMessage(L"Display Ux");
        m_current_max = m_U_max;
        return m_Ux;
        break;
    case 1:
        m_pDraw->setMessage(L"Display Uy");
        m_current_max = m_U_max;
        return m_Uy;
        break;
    case 3:
        m_pDraw->setMessage(L"Dislapy P");
        m_current_max = m_p_max;
        return m_p;
        break;
    default:
        m_pDraw->setMessage(L"Done");
    }

    return nullptr;
}
void Test::drawDisplayFrames() {
    double* data_to_draw = getFrameToDisplay();
    if(data_to_draw!=nullptr)
        m_pDraw->drawData(data_to_draw, m_current_max);
}
GenImage* Test::handleMouse() {
    drawDisplayFrames();
    return m_pCurrent_GenImage;
}
void Test::find_max(double& max, const double* data) {
    max = 0.0;
    int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    for (int i = 0; i < size; i++) {
        double abs_dat = std::abs(data[i]);
        if (abs_dat > max)
            max = abs_dat;
    }
}
void Test::find_max() {
    find_max(m_Ux_max, m_Ux);
    find_max(m_Uy_max, m_Uy);
    m_U_max = (m_Ux_max > m_Uy_max) ? m_Ux_max : m_Uy_max;
    find_max(m_p_max, m_p);
    m_current_max = 0.0;
}

void Test::send_frame_to_host(double* pFrame, const double* dev_data) {
    int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    cudaError_t cudaStatus = cudaMemcpy(pFrame, dev_data, size * sizeof(double), cudaMemcpyDeviceToHost);
}
void Test::send_frame_to_host(double* pFrame, const double* dev_data, int dat_len) {
    cudaError_t cudaStatus = cudaMemcpy(pFrame, dev_data, dat_len * sizeof(double), cudaMemcpyDeviceToHost);
}