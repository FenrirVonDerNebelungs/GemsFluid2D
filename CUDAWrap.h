#pragma once
#ifndef CUDAWRAP_H
#define CUDAWRAP_H
#ifndef BASE_H
#include "Base.h"
#endif

class Test;

class CUDAWrap {
public:
    CUDAWrap(
        int blocks_side_dim = 4,//8,
        int threads_side_dim = 16,
        double in_delta_t = 1e-3,
        double in_delta_x = 1e-3,
        double in_nu = 1.0,
        int jacobi_minBlocks_side_dim=2, /*must be such that 2^some power * minBlocks_side dim = blocks_side_dim */
		int jacobi_minThreads_side_dim = 4, /*must be such that 2^some power * minThreads_side dim = threads_side_dim */
        int in_max_jacobi_loops = 200
    );
    ~CUDAWrap();

    int runCUDA(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames, double jacobi_filter[]);

	s_WH getGridWidthHeight() { s_WH wh; wh.width = grid_width; wh.height = grid_height; return wh; }

	friend class Test;
protected:

    double* m_dev_Ux[2];/* 2 frames for in and out, each frame is grid_width*grid_height doubles */
    double* m_dev_Uy[2];/* 2 frames for in and out, each frame is grid_width*grid_height doubles */
    double* m_dev_p[2];/* 2 frames for in and out, each frame is grid_width*grid_height doubles */
    double* m_dev_scratch;

    double* m_filter;

    int numBlocks_side;
    int numThreads_side;
    int numBlocks_for_1D;
	int numThreads_for_1D;
    int grid_width;
    int grid_height;
    double delta_t;
    double delta_x;
    double nu;

    int max_jacobi_loops;
    int jacobi_stack_height;
    int jacobi_block_expansion_pow;
	int jacobi_thread_expansion_pow;
	int jacobi_minBlocks_side;
	int jacobi_minThreads_side;
	int* jacobi_scratch_stack_sizes;
    int* jacobi_stack_numBlocks;
	int* jacobi_stack_numThreads;
    s_WH* jacobi_stack_WH;
    double** jacobi_scratch_stack;
    double* jacobi_scratch;
    double** b_stack;
    double** Wx_stack;
    double** Wy_stack;
    double* jacobi_alpha;
	double  jacobi_rbeta;
    double* jacobi_delta_x;

    int runNV(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames, double jacobi_filter[]);

	void runFrame(double* Ux[], double* Uy[], double* p[], double* scratch, int& frame_index, int& p_frame_index, s_force& force);

    void advection(double* Ux[], double* Uy[], int frame_index);
    void advection_backtrace(double* relPos_i, double* relPos_j, /*const*/ double* Ux[], /*const*/ double* Uy[], int frame_index);/* test function */
    void viscous_diffusion(double* Ux[], double* Uy[], int frame_index);
    void apply_force(double* Ux[], double* Uy[], int frame_index, s_force& force);
    void compute_pressure(double* Ux[], double* Uy[], double* p[], double* scratch, int frame_index, int p_frame_index);
    void subtract_pressure_gradient(double* Ux[], double* Uy[], double* p[], int frame_index, int p_frame_index);

    void jacobi_fill_stacks(
        const double* frame_in, 
        const double* b,
        const double* Wx,
		const double* Wy);
    void jacobi_send_frame_down_stack(
        const double* hi_frame, 
        int hi_frame_width_height, 
        int lo_stack_index);

    void jacobi_frame(
        double* frame_out, 
        const double* frame_in, 
        const double* b, 
        const double* Wx, 
        const double* Wy, 
        const double& alpha, 
		const double& rbeta,
		const int& numBlocks_s,
		const int& numThreads_s);

    void jacobi_loop(
        double* X[], /*frame in is X[0], frame out is also X[0], X[1] in is swap, the pointer X[0] may = the oringinal pointer X[1]*/
        const double* b, 
        double alpha, 
        double rbeta, 
        int numBlocks_s,
		int numThreads_s,
        const double* Wx=nullptr, 
        const double* Wy=nullptr);
    void jacobi_run(
		double* X[],
        const double* b,
        const double* Wx,
        const double* Wy,
        int frame_index);
	void divergence(double* div, const double* Ux, const double* Uy);

	void gradient(double* dp_dx, double* dp_dy, double* p[], int p_frame_index);
    void laplacian(double* lap, double* p[], double* scratch_x, double* scratch_y, int p_frame_index);
	void bilinearAprox_scaledFrame(double* Ux_scaled, double* Uy_scaled, const double* Ux, const double* Uy, int scale_factor=6);
	void copyMemory_for_standard_grid(double* dest, const double* src);

    /* init/release  */
    bool initDevMem(int size, double jacobi_filter[]);
	void releaseDevMem();
    /*utility*/
    int findLog_base2(int val);
    int find2Pow(int pow);
	int find_stack_WH_and_redFactor(s_WH& wh, int stack_level, int stack_height);
	void find_stack_BlocksNThreads(int& numBlocks_s, int& numThreads_s, int reduction_factor);
};

#endif