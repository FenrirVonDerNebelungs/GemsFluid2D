#include "Test.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
Test::Test(int num_display_frames) :m_Ux_max(0.0), m_Uy_max(0.0), m_p_max(0.0), m_current_max(0.0), m_mouse_clicks(0), m_num_display_frames(num_display_frames), m_max_error(0.0) {
    m_pCUDA_wrap = new CUDAWrap();
	int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
	m_Ux = new double[size];
	m_Uy = new double[size];
	m_p = new double[size];
	m_scratch = new double[size];
	m_display_frames = new double* [m_num_display_frames];
    for (int i = 0; i < m_num_display_frames; i++) {
		m_display_frames[i] = new double[size];
    }
	m_pFluid_animate = new FluidAnimate();
	s_WH wh = getGridWidthHeight();
	m_pFluid_animate->init(wh);
	m_force = m_pFluid_animate->getForce(0);
}
Test::~Test() {
	if (m_pFluid_animate != nullptr)
		delete m_pFluid_animate;
    if (m_display_frames != nullptr) {
        for (int i = 0; i < m_num_display_frames; i++) {
			if (m_display_frames[i] != nullptr)
                delete[] m_display_frames[i];
        }
		delete[] m_display_frames;
    }
	if (m_scratch != nullptr)
		delete[] m_scratch;
    if(m_p!= nullptr)
		delete[] m_p;
    if(m_Uy!=nullptr)
		delete[] m_Uy;
	if (m_Ux != nullptr)
		delete[] m_Ux;
	if (m_pCUDA_wrap != nullptr)
        delete m_pCUDA_wrap;
}
void Test::runTestFrame(double* Ux[], double* Uy[], double* p[], double* scratch, int& frame_index, int& p_frame_index, s_force& force) {
    int frame_in_index = frame_index;
    int p_frame_in_index = p_frame_index;
    if (force.active) {
        m_pCUDA_wrap->apply_force(Ux, Uy, frame_in_index, force);
        reverseFrameIndex(frame_in_index);
    }
    fill_display_frame(m_Ux, Ux[frame_in_index]);
	fill_display_frame(m_Uy, Uy[frame_in_index]);
	m_pCUDA_wrap->divergence(scratch, Ux[frame_in_index], Uy[frame_in_index]);
	fill_display_frame(m_scratch, scratch);
    static double alpha = -m_pCUDA_wrap->delta_x * m_pCUDA_wrap->delta_x;
    static double rbeta = 0.25;
	runJacobiTest(p, scratch, p_frame_in_index, alpha, rbeta);
	reverseFrameIndex(p_frame_in_index);
	fill_display_frame(m_p, p[p_frame_in_index]);
	frame_index = frame_in_index;
	p_frame_index = p_frame_in_index;
}
void Test::runJacobiTest(double* U[], double* b, int frame_index, double alpha, double rbeta) {
    dim3 numBlocks(m_pCUDA_wrap->numBlocks_side, m_pCUDA_wrap->numBlocks_side);
    dim3 numThreads(m_pCUDA_wrap->numThreads_side, m_pCUDA_wrap->numThreads_side);
    cudaError_t cudaStatus = cudaSuccess;
    int original_frame_index = frame_index;
    s_frame_index frame_i = getFrameIndex(frame_index);
    int num_jacobi_loops = 0;
    do {
		m_pCUDA_wrap->jacobi_frame(U[frame_i.out], U[frame_i.in], b, alpha, rbeta);
		fill_display_frame(num_jacobi_loops, U[frame_i.out]);
        swapFrameIndexes(frame_i);
        num_jacobi_loops++;
    } while (num_jacobi_loops < m_pCUDA_wrap->max_jacobi_loops && num_jacobi_loops<m_num_display_frames);
    fixFramePointers(U, frame_i, original_frame_index);/* set Ux so that frame_out_index will point to the Ux results of the jacobi loop*/
}
int Test::runTest(int sim_frames){
    if (sim_frames > m_num_display_frames)
        return 1;
	int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    std::memset(m_Ux, 0, size * sizeof(double));
    std::memset(m_Uy, 0, size * sizeof(double));
    std::memset(m_p, 0, size * sizeof(double));
    return runCUDA(m_Ux, m_Uy, m_p, m_force, sim_frames);
}
int Test::runCUDA(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames) {
    unsigned int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
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
            runTestFrame(dev_Ux, dev_Ux, dev_p, scratch, frame_index, p_frame_index, force);
            frames_run++;
        } while (frames_run < sim_frames);
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess)
            fprintf(stderr, "Cuda kernel launches failed:%s\n", cudaGetErrorString(cudaStatus));
    }
    /*
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(Ux, dev_Ux[0], size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(Uy, dev_Uy[0], size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(pressure, dev_p[0], size * sizeof(double), cudaMemcpyDeviceToHost);
    */
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
double* Test::getFrameToDisplay() {
    int mouse_clicks = m_mouse_clicks;
    m_mouse_clicks++;
    switch (mouse_clicks) {
    case 0:
        return m_Ux;
        break;
    case 1:
        return m_Uy;
        break;
    case 2:
        return m_scratch;
        break;
    case 3:
        find_max();
    }
	if (mouse_clicks >= 3 && mouse_clicks < 3 + m_num_display_frames)
    {
        int extra_clicks = mouse_clicks - 3;
        if (extra_clicks < m_num_display_frames)
            return m_display_frames[extra_clicks];
	}
    else if (mouse_clicks < 3 + 2 * m_num_display_frames) {
        int extra_clicks = mouse_clicks - 3 - m_num_display_frames;
		if (extra_clicks == 0)
			set_display_frames_to_errors();
        if (extra_clicks < m_num_display_frames) {
            return m_display_frames[extra_clicks];
        }
    }
    
    return m_p;
}
void Test::find_max(double& max, const double* data) {
    max = 0.0;
    int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    for (int i = 0; i < size; i++) {
        if (data[i] > max)
            max = data[i];
    }
}
void Test::find_max() {
    m_current_max = 0.0;
    for(int i=0; i<m_num_display_frames; i++) {
        double frame_max = 0.0;
        find_max(frame_max, m_display_frames[i]);
        if(frame_max > m_current_max)
			m_current_max = frame_max;
	}
}
double Test::findErrors() {
    set_display_frames_to_errors();
    return m_max_error;
}
void Test::set_display_frames_to_errors() {
    int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
	m_max_error = 0.0;
    for (int i = 1; i < m_num_display_frames; i++) {
		double display_frame_max_error = 0.0;
        for (int j = 0; j < size; j++) {
            m_display_frames[(i-1)][j] = m_display_frames[i][j] - m_display_frames[(i-1)][j];
            if(m_display_frames[(i-1)][j]>display_frame_max_error)
				display_frame_max_error = m_display_frames[(i - 1)][j];
        }
		if (display_frame_max_error > m_max_error)
			m_max_error = display_frame_max_error;
    }
}
void Test::fill_display_frame(int frame_index, const double* dev_data) {
    int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    cudaError_t cudaStatus = cudaMemcpy(m_display_frames[frame_index], dev_data, size * sizeof(double), cudaMemcpyDeviceToHost);
}
void Test::fill_display_frame(double* pFrame, const double* dev_data) {
    int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    cudaError_t cudaStatus = cudaMemcpy(pFrame, dev_data, size * sizeof(double), cudaMemcpyDeviceToHost);
}