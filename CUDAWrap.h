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
        int in_max_jacobi_loops = 400
    );
    ~CUDAWrap();

    int runCUDA(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames);

	s_WH getGridWidthHeight() { s_WH wh; wh.width = grid_width; wh.height = grid_height; return wh; }

	friend class Test;
protected:
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
    double** jacobi_scratch_stack;
    double* jacobi_scratch;
    double** b_stack;

    int runNV(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames);

	void runFrame(double* Ux[], double* Uy[], double* p[], double* scratch, int& frame_index, int& p_frame_index, s_force& force);

    void advection(double* Ux[], double* Uy[], int frame_index);
    void advection_backtrace(double* relPos_i, double* relPos_j, /*const*/ double* Ux[], /*const*/ double* Uy[], int frame_index);/* test function */
    void viscous_diffusion(double* Ux[], double* Uy[], double* scratch, int frame_index);
    void apply_force(double* Ux[], double* Uy[], int frame_index, s_force& force);
    void compute_pressure(double* Ux[], double* Uy[], double* p[], double* scratch, int frame_index, int p_frame_index);
    void subtract_pressure_gradient(double* Ux[], double* Uy[], double* p[], int frame_index, int p_frame_index);

    void jacobi_b_stack(const double* b);
    void jacobi_run_stack(
        const double* frame_in,
        const double& alpha,
        const double& rbeta);
    void jacobi_base(
        double* frame_out,
        const double* frame_in,
        const double* b,
        const double* Wx,
        const double* Wy,
        const double& alpha,
        const double& rbeta);

    void jacobi_frame(
        double* frame_out, 
        const double* frame_in, 
        const double* b, 
        const double* Wx, 
        const double* Wy, 
        const double& alpha, 
		const double& rbeta);



    void jacobi_loop(
        double* X[], 
        double* b, 
        int frame_index,
        double alpha, 
        double rbeta, 
        const double* Wx=nullptr, const double* Wy=nullptr);
	void divergence(double* div, const double* Ux, const double* Uy);

	void bilinearAprox_scaledFrame(double* Ux_scaled, double* Uy_scaled, const double* Ux, const double* Uy, int scale_factor=6);

    /*utility*/
    int findPow2(int val);
	bool find_reduced_BlocksNThreads(int& numBlocks_s, int& numThreads_s, int red_factor);
};

#endif