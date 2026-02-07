#include "Test.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
Test::Test(int num_display_frames) :m_num_display_frames(num_display_frames) {
    m_pCUDA_wrap = new CUDAWrap();
	int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
	m_Ux = new double[size];
	m_Uy = new double[size];
	m_p = new double[size];
	m_display_frames = new double* [m_num_display_frames];
    for (int i = 0; i < m_num_display_frames; i++) {
		m_display_frames[i] = new double[size];
    }
}
Test::~Test() {
    if (m_display_frames != nullptr) {
        for (int i = 0; i < m_num_display_frames; i++) {
			if (m_display_frames[i] != nullptr)
                delete[] m_display_frames[i];
        }
		delete[] m_display_frames;
    }
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
    dim3 numBlocks(m_pCUDA_wrap->numBlocks_side, m_pCUDA_wrap->numBlocks_side);
    dim3 numThreads(m_pCUDA_wrap->numThreads_side, m_pCUDA_wrap->numThreads_side);
    /* jacobi with \alpha = -\deltax^2 and b = \frac{1}{\deltat} \nabla \cdot \vec{u} and \beta = 4 */
    static double alpha = -m_pCUDA_wrap->delta_x * m_pCUDA_wrap->delta_x;
    static double rbeta = 0.25;
    static double inv_2delta_x = 1.0 / (2.0 * m_pCUDA_wrap->delta_x);
    divergence << <numBlocks, numThreads >> > (scratch, Ux[frame_index], Uy[frame_index], inv_2delta_x, grid_width, grid_height);
    frame_index = frame_in_index;
    p_frame_index = p_frame_in_index;
}
void Test::runJacobiTest(double* U[], double* b, int frame_index, double alpha, double rbeta) {

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