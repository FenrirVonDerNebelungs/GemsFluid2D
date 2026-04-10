#include "Test.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

Test::Test(int num_display_frames, int blow_factor) :
    m_current_frame(0),
    m_blow_factor(blow_factor),
    m_Ux_max(0.0),
    m_Uy_max(0.0),
    m_U_max(0.0),
    m_div_max(0.0),
    m_Ux_new_max(0.0),
    m_Uy_new_max(0.0),
    m_U_new_max(0.0),
    m_p_max(0.0),
    m_current_max(0.0),
    m_mouse_clicks(0),
    m_num_display_frames(num_display_frames),
    m_max_error(0.0),
    m_image_sup(0.0)
{
    m_pCUDA_wrap = new CUDAWrap();
	int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
	m_Ux = new double[size];
	m_Uy = new double[size];
    m_relPos_x = new double[size];
    m_relPos_y = new double[size];
	m_p = new double[size];
	m_scratch = new double[size];
	m_display_frames = new double* [m_num_display_frames];
    for (int i = 0; i < m_num_display_frames; i++) {
		m_display_frames[i] = new double[size];
    }
    m_Ux_new = new double[size];
    m_Uy_new = new double[size];
    m_U_div = new double[size];
    int blown_size = size * m_blow_factor * m_blow_factor;
    m_Ux_bilinear = new double[blown_size];
    m_Uy_bilinear = new double[blown_size];
    m_pPyTrans = new PyTrans();
	m_pFluid_animate = new FluidAnimate();
	s_WH wh = getGridWidthHeight();
	m_pFluid_animate->init(wh);
	m_force = m_pFluid_animate->getForce(0);
    m_drawVelocity = new drawVelocity(m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_pCUDA_wrap->delta_x, m_pCUDA_wrap->delta_t, nullptr, m_blow_factor);
    m_drawTest = new drawTest(m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height);
    m_pCurrent_GenImage = m_drawTest->m_pGenImage;
}
Test::~Test() {
    m_pCurrent_GenImage = nullptr;
    if (m_drawTest != nullptr)
        delete m_drawTest;
    if (m_drawVelocity != nullptr)
        delete m_drawVelocity;
    if (m_Uy_bilinear != nullptr)
        delete [] m_Uy_bilinear;
    if (m_Ux_bilinear != nullptr)
        delete[] m_Ux_bilinear;
    if (m_U_div != nullptr)
        delete[] m_U_div;
    if (m_Uy_new != nullptr)
        delete[] m_Uy_new;
    if (m_Ux_new != nullptr)
        delete[] m_Ux_new;
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
    if (m_relPos_x != nullptr)
        delete[] m_relPos_x;
    if (m_relPos_y != nullptr)
        delete[] m_relPos_y;
    if(m_Uy!=nullptr)
		delete[] m_Uy;
	if (m_Ux != nullptr)
		delete[] m_Ux;
	if (m_pCUDA_wrap != nullptr)
        delete m_pCUDA_wrap;
    if (m_pPyTrans != nullptr)
        delete m_pPyTrans;
}
void Test::runTestToPressure(double* Ux[], double* Uy[], double* p[], double* scratch, int& frame_index, int& p_frame_index, s_force& force) {
    int frame_in_index = frame_index;
    int p_frame_in_index = p_frame_index;
    if (force.active) {
        m_pCUDA_wrap->apply_force(Ux, Uy, frame_in_index, force);
        reverseFrameIndex(frame_in_index);
    }
    fill_display_frame(m_Ux, Ux[frame_in_index]);
	fill_display_frame(m_Uy, Uy[frame_in_index]);
    m_pPyTrans->cacheGrid(m_Ux, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::W_code, n_PyTrans::X_code, n_PyTrans::after_force_code);
    m_pPyTrans->cacheGrid(m_Uy, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::W_code, n_PyTrans::Y_code, n_PyTrans::after_force_code);
	m_pCUDA_wrap->divergence(scratch, Ux[frame_in_index], Uy[frame_in_index]);
	fill_display_frame(m_scratch, scratch);
    m_pPyTrans->cacheGrid(m_scratch, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::DivW_code);
    static double alpha = -m_pCUDA_wrap->delta_x * m_pCUDA_wrap->delta_x;
    static double rbeta = 0.25;
	runJacobiTest(p, scratch, p_frame_in_index, alpha, rbeta, Ux[frame_in_index], Uy[frame_in_index]);
	reverseFrameIndex(p_frame_in_index);
	fill_display_frame(m_p, p[p_frame_in_index]);
    m_pPyTrans->cacheGrid(m_p, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::P_code, n_PyTrans::Scalar_code, n_PyTrans::end_frame_code);
    m_pCUDA_wrap->subtract_pressure_gradient(Ux, Uy, p, frame_in_index, p_frame_in_index);
    reverseFrameIndex(frame_in_index);
    fill_display_frame(m_Ux_new, Ux[frame_in_index]);
    fill_display_frame(m_Uy_new, Uy[frame_in_index]);
    m_pPyTrans->cacheGrid(m_Ux, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::U_code, n_PyTrans::X_code, n_PyTrans::end_frame_code);
    m_pPyTrans->cacheGrid(m_Uy, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::U_code, n_PyTrans::Y_code, n_PyTrans::end_frame_code);
    m_pCUDA_wrap->divergence(scratch, Ux[frame_in_index], Uy[frame_in_index]);//this should be zero
    fill_display_frame(m_U_div, scratch);
    m_pPyTrans->cacheGrid(m_U_div, m_pCUDA_wrap->grid_width, m_pCUDA_wrap->grid_height, m_current_frame, n_PyTrans::DivU_code, n_PyTrans::Scalar_code, n_PyTrans::end_frame_code);
	frame_index = frame_in_index;
	p_frame_index = p_frame_in_index;
}
void Test::runTestAdvection(
    double* Ux[], 
    double* Uy[], 
    double* p[], 
    double* scratch, 
    double* Ux_bilinear, 
    double* Uy_bilinear, 
    double* relPos_x,
    double* relPos_y,
    int& frame_index, 
    int& p_frame_index, 
    s_force& force) 
{
    cudaError_t cudaStatus = cudaSuccess;
    int frame_in_index = frame_index;
    int p_frame_in_index = p_frame_index;
    runTestToPressure(Ux, Uy, p, scratch, frame_in_index, p_frame_in_index, force); /*end stuff is in in-indexes*/
    p_frame_index = p_frame_in_index;
    fill_display_frame(m_Ux, Ux[frame_in_index]);
    fill_display_frame(m_Uy, Uy[frame_in_index]);
    size_t size_in_pix_of_blown_image = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height * m_drawVelocity->getBlowFactor() * m_drawVelocity->getBlowFactor();
    size_t size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;

    m_pCUDA_wrap->bilinearAprox_scaledFrame(Ux_bilinear, Uy_bilinear, Ux[frame_in_index], Uy[frame_in_index], m_drawVelocity->getBlowFactor());
    fill_display_frame(m_Ux_bilinear, Ux_bilinear, static_cast<int>(size_in_pix_of_blown_image));
    fill_display_frame(m_Uy_bilinear, Uy_bilinear, static_cast<int>(size_in_pix_of_blown_image));
    /*m_pPyTrans->cacheGrid(
        m_Ux_bilinear, 
        (m_pCUDA_wrap->grid_width*m_drawVelocity->getBlowFactor()), 
        (m_pCUDA_wrap->grid_height*m_drawVelocity->getBlowFactor()), 
        m_current_frame, 
        n_PyTrans::U_code, 
        n_PyTrans::X_code, 
        n_PyTrans::start_frame_code, 
        0, 
        m_drawVelocity->getBlowFactor()
    );
    m_pPyTrans->cacheGrid(
        m_Uy_bilinear,
        (m_pCUDA_wrap->grid_width * m_drawVelocity->getBlowFactor()),
        (m_pCUDA_wrap->grid_height * m_drawVelocity->getBlowFactor()),
        m_current_frame,
        n_PyTrans::U_code,
        n_PyTrans::Y_code,
        n_PyTrans::start_frame_code,
        0,
        m_drawVelocity->getBlowFactor()
    );*/
    m_pCUDA_wrap->advection_backtrace(relPos_x, relPos_y, Ux, Uy, frame_in_index);
    fill_display_frame(m_relPos_x, relPos_x);
    fill_display_frame(m_relPos_y, relPos_y);
    m_pPyTrans->cacheGrid(
        m_relPos_x,
        (m_pCUDA_wrap->grid_width),
        (m_pCUDA_wrap->grid_height),
        m_current_frame,
        n_PyTrans::relPos_code,
        n_PyTrans::X_code
    );
    m_pPyTrans->cacheGrid(
        m_relPos_y,
        (m_pCUDA_wrap->grid_width),
        (m_pCUDA_wrap->grid_height),
        m_current_frame,
        n_PyTrans::relPos_code,
        n_PyTrans::Y_code
    );
    m_pCUDA_wrap->advection(Ux, Uy, frame_in_index);
    s_frame_index s_frame = getFrameIndex(frame_in_index);
    fill_display_frame(m_Ux_new, Ux[s_frame.out]);
    fill_display_frame(m_Uy_new, Uy[s_frame.out]);
    m_pPyTrans->cacheGrid(
        m_Ux_new,
        (m_pCUDA_wrap->grid_width),
        (m_pCUDA_wrap->grid_height),
        m_current_frame,
        n_PyTrans::U_code,
        n_PyTrans::X_code, 
        n_PyTrans::after_advection_code
    );
    m_pPyTrans->cacheGrid(
        m_Uy_new,
        (m_pCUDA_wrap->grid_width),
        (m_pCUDA_wrap->grid_height),
        m_current_frame,
        n_PyTrans::U_code,
        n_PyTrans::Y_code,
        n_PyTrans::after_advection_code
    );
    frame_index = s_frame.out;
}
void Test::runTestFrame(
    double* Ux[], 
    double* Uy[], 
    double* p[], 
    double* scratch, 
    double* Ux_bilinear, 
    double* Uy_bilinear, 
    double* relPos_x,
    double* relPos_y,
    int& frame_index, 
    int& p_frame_index, 
    s_force& force) {
    //runTestToPressure(Ux, Uy, p, scratch, frame_index, p_frame_index, force);
    runTestAdvection(Ux, Uy, p, scratch, Ux_bilinear, Uy_bilinear, relPos_x, relPos_y, frame_index, p_frame_index, force);
}
void Test::runJacobiTest(double* U[], double* b, int frame_index, double alpha, double rbeta, const double* Wx, const double* Wy) {
    dim3 numBlocks(m_pCUDA_wrap->numBlocks_side, m_pCUDA_wrap->numBlocks_side);
    dim3 numThreads(m_pCUDA_wrap->numThreads_side, m_pCUDA_wrap->numThreads_side);
    cudaError_t cudaStatus = cudaSuccess;
    int original_frame_index = frame_index;
    s_frame_index frame_i = getFrameIndex(frame_index);
    int num_jacobi_loops = 0;
    int num_display_frames = 0;
    do {
		m_pCUDA_wrap->jacobi_frame(U[frame_i.out], U[frame_i.in], b, Wx, Wy, alpha, rbeta);
        if ( ((num_jacobi_loops+1) % 10) == 0 && num_display_frames<m_num_display_frames) {
            fill_display_frame(num_display_frames, U[frame_i.out]);
            num_display_frames++;
        }
        swapFrameIndexes(frame_i);
        num_jacobi_loops++;
    } while (num_jacobi_loops < m_pCUDA_wrap->max_jacobi_loops);
    fixFramePointers(U, frame_i, original_frame_index);/* set Ux so that frame_out_index will point to the Ux results of the jacobi loop*/
}

int Test::runTest(int sim_frames){
    if (sim_frames > m_num_display_frames)
        return 1;
	int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    std::memset(m_Ux, 0, size * sizeof(double));
    std::memset(m_Uy, 0, size * sizeof(double));
    std::memset(m_p, 0, size * sizeof(double));
    runCUDA(m_Ux, m_Uy, m_p, m_force, sim_frames);
    m_drawVelocity->init(m_Ux, m_Uy, m_relPos_x, m_relPos_y, m_Ux_new, m_Uy_new, m_Ux_bilinear, m_Uy_bilinear);
    find_max();
    return 0;
}
int Test::runCUDA(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames) {
    unsigned int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    unsigned int size_in_pix_of_blown_image = size * m_blow_factor * m_blow_factor;
    double* dev_Ux[] = { 0,0 };/*two frames */
    double* dev_Uy[] = { 0,0 };
    double* dev_p[] = { 0,0 };
    double* scratch = 0;
    double* dev_Ux_bilinear = 0;
    double* dev_Uy_bilinear = 0;
    double* dev_relPos_x = 0; 
    double* dev_relPos_y = 0;

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
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&dev_Ux_bilinear, size_in_pix_of_blown_image * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&dev_Uy_bilinear, size_in_pix_of_blown_image * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&dev_relPos_x, size * sizeof(double));
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMalloc((void**)&dev_relPos_y, size * sizeof(double));
    if (cudaStatus != cudaSuccess)
        fprintf(stderr, "cudaMalloc failed!");


    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(dev_Ux[0], Ux, size * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus == cudaSuccess)
        cudaStatus = cudaMemcpy(dev_Uy[0], Uy, size * sizeof(double), cudaMemcpyHostToDevice);

    int grid_stream_len = m_pCUDA_wrap->getGridWidthHeight().height * m_pCUDA_wrap->getGridWidthHeight().width;
    //int grid_exp_stream_len = grid_stream_len * m_blow_factor * m_blow_factor;
    int total_grid_stream_len = grid_stream_len * (11 + m_num_display_frames);// +grid_exp_stream_len * 2;
    int total_num_stream_headers = (11 + m_num_display_frames); //+ 2;
    m_pPyTrans->init("Dat/frames.dat", total_grid_stream_len, total_num_stream_headers);
    if (cudaStatus == cudaSuccess) {
        int frames_run = 0;
        int frame_index = 0;
        int p_frame_index = 0;
        do {
            m_current_frame = frames_run;
            runTestFrame(
                dev_Ux, 
                dev_Ux, 
                dev_p, 
                scratch, 
                dev_Ux_bilinear, 
                dev_Uy_bilinear, 
                dev_relPos_x, 
                dev_relPos_y, 
                frame_index, 
                p_frame_index, 
                force);
            m_pPyTrans->resetAndWrite();
            if (frames_run > 2)
                force.active = false;
            frames_run++;
        } while (frames_run < sim_frames);
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess)
            fprintf(stderr, "Cuda kernel launches failed:%s\n", cudaGetErrorString(cudaStatus));
    }
    m_pPyTrans->releaseAndWrite();
    /*debug*/
    //m_pPyTrans->init("Dat/test.dat", total_grid_stream_len, total_num_stream_headers);
    //double test_py_out[4] = { 0.32, -1.45, 12.34e-10, -59.0e7 };
    //m_pPyTrans->cacheDStream(test_py_out, 4);
    //m_pPyTrans->releaseAndWrite();
    /*     */
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

    cudaFree(dev_relPos_x);
    cudaFree(dev_relPos_y);
    cudaFree(dev_Uy_bilinear);
    cudaFree(dev_Ux_bilinear);
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
        m_pCurrent_GenImage = m_drawTest->getGenImage();
        m_drawTest->setMessage(L"Display Ux after force");
        m_image_sup = m_U_max;
        return m_Ux;
        break;
    case 1:
        m_drawTest->setMessage(L"Display Uy after force");
        m_image_sup = m_U_max;
        return m_Uy;
        break;
    case 2:
        m_drawTest->setMessage(L"Divergence of U after force");
        m_image_sup = m_div_max;
        return m_scratch;
        break;
    case 3:
        m_image_sup = m_current_max;
        m_drawTest->setMessage(L"Displaying jacobi convergence to delta^2 p = div U");
    }
    int num_init_cases = 3;
	if (mouse_clicks >= num_init_cases && mouse_clicks < num_init_cases + m_num_display_frames)
    {
        int extra_clicks = mouse_clicks - num_init_cases;
        if (extra_clicks >= 1)
            m_drawTest->setMessage(L"");
        if (extra_clicks < m_num_display_frames)
            return m_display_frames[extra_clicks];
	}
    else if (mouse_clicks < num_init_cases + 2 * m_num_display_frames) {
        int extra_clicks = mouse_clicks - num_init_cases - m_num_display_frames;
        if (extra_clicks == 0) {
            set_display_frames_to_errors();
            m_image_sup = m_max_error;
            m_drawTest->setMessage(L"displaying errors between jacobi frames");
        }
        else
            m_drawTest->setMessage(L"");
        if (extra_clicks < (m_num_display_frames-1)) {
            return m_display_frames[extra_clicks];
        }
        else {
            m_drawTest->setMessage(L"last error frame");
            return m_display_frames[m_num_display_frames - 2];
        }
    }
    int second_batch_mouse_clicks = mouse_clicks - (num_init_cases + 2 * m_num_display_frames);
    switch (second_batch_mouse_clicks) {
    case 0:
        m_drawTest->setMessage(L"displaying computed pressure");
        m_image_sup = m_p_max;
        return m_p;
        break;
    case 1:
        m_drawTest->setMessage(L"new velocity Ux");
        m_image_sup = m_U_new_max;
        return m_Ux_new;
        break;
    case 2:
        m_drawTest->setMessage(L"new velocity Uy");
        m_image_sup = m_U_new_max;
        return m_Uy_new;
        break;
    case 3:
        m_drawTest->setMessage(L"divergence of U after subtraction of pressure gradient");
        m_image_sup = m_U_new_max;
        return m_U_div;
    }
    int total_second_batch_mouse_clicks = 4;
    int advection_mouse_clicks = mouse_clicks - ((num_init_cases+2*m_num_display_frames) + total_second_batch_mouse_clicks);
    switch (advection_mouse_clicks) {
    case 0:
        m_pCurrent_GenImage = m_drawVelocity->getGenImage();
        m_drawTest->setMessage(L"Divergenceless velocity after pressure grad subtracted");
        m_drawVelocity->setImage_to_velocity();
        return nullptr;
    case 1:
        m_drawTest->setMessage(L"Velocity drawn with arrow shifted backwards");
        m_drawVelocity->setImage_to_backVelocity();
        return nullptr;
    case 2:
        m_drawTest->setMessage(L"Back displacement from velocity estimated by advection code");
        m_drawVelocity->setImage_to_displacement();
        return nullptr;
    case 3:
        m_drawTest->setMessage(L"New velocity estimated from bilinear approx by code");
        m_drawVelocity->setImage_to_newVelocity();
        return nullptr;
    }
    if (advection_mouse_clicks > 3 && advection_mouse_clicks <= 6) {
        int index_of_bi = advection_mouse_clicks - 4;
        m_drawTest->setMessage(L"Bilinear Ux");
        m_drawVelocity->setImage_to_bilinear_Ux(index_of_bi);
        return nullptr;
    }
    else if (advection_mouse_clicks <= 9) {
        int index_of_bi = advection_mouse_clicks - 7;
        m_drawTest->setMessage(L"Bilinear Uy");
        m_drawVelocity->setImage_to_bilinear_Uy(index_of_bi);
        return nullptr;
    }
    m_drawTest->setMessage(L"Done");
    return nullptr;
}
void Test::drawDisplayFrames() {
    double* data_to_draw = getFrameToDisplay();
    if(data_to_draw!=nullptr)
        m_drawTest->drawData(data_to_draw, m_image_sup);
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
    find_max(m_div_max, m_scratch);
    find_max(m_Ux_new_max, m_Ux_new);
    find_max(m_Uy_new_max, m_Uy_new);
    m_U_new_max = (m_Ux_new_max > m_Uy_new_max) ? m_Ux_new_max : m_Uy_new_max;
    /*debug*/
    find_max(m_U_new_max, m_U_div);
    /*     */
    find_max(m_p_max, m_p);
    m_current_max = 0.0;
    for(int i=0; i<m_num_display_frames; i++) {
        double frame_max = 0.0;
        find_max(frame_max, m_display_frames[i]);
        if(frame_max > m_current_max)
			m_current_max = frame_max;
	}
}

void Test::set_display_frames_to_errors() {
    int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
	m_max_error = 0.0;
    for (int i = 1; i < m_num_display_frames; i++) {
		double display_frame_max_error = 0.0;
        for (int j = 0; j < size; j++) {
            m_display_frames[(i-1)][j] = m_display_frames[i][j] - m_display_frames[(i-1)][j];
            double abs_error = abs(m_display_frames[(i - 1)][j]);
            if(abs_error>display_frame_max_error)
				display_frame_max_error = abs_error;
        }
		if (display_frame_max_error > m_max_error)
			m_max_error = display_frame_max_error;
    }
}
void Test::fill_display_frame(int frame_index, const double* dev_data) {
    int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    cudaError_t cudaStatus = cudaMemcpy(m_display_frames[frame_index], dev_data, size * sizeof(double), cudaMemcpyDeviceToHost);
    m_pPyTrans->cacheGrid(
        m_display_frames[frame_index],
        (m_pCUDA_wrap->grid_width),
        (m_pCUDA_wrap->grid_height),
        m_current_frame,
        n_PyTrans::P_code,
        n_PyTrans::Scalar_code,
        n_PyTrans::after_force_code,
        frame_index
    );
}
void Test::fill_display_frame(double* pFrame, const double* dev_data) {
    int size = m_pCUDA_wrap->grid_width * m_pCUDA_wrap->grid_height;
    cudaError_t cudaStatus = cudaMemcpy(pFrame, dev_data, size * sizeof(double), cudaMemcpyDeviceToHost);
}
void Test::fill_display_frame(double* pFrame, const double* dev_data, int dat_len) {
    cudaError_t cudaStatus = cudaMemcpy(pFrame, dev_data, dat_len * sizeof(double), cudaMemcpyDeviceToHost);
}