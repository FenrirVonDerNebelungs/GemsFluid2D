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
        int blocks_side_dim = 8,
        int threads_side_dim = 16,
        double delta_t = 1e-3,
        double delta_x = 1e-3,
        double nu = 1.0,
        int max_jacobi_loops = 3
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

    int runNV(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames);

	void runFrame(double* Ux[], double* Uy[], double* p[], double* scratch, int& frame_index, int& p_frame_index, s_force& force);
	void runTestFrame(double* Ux[], double* Uy[], double* p[], double* scratch, int& frame_index, int& p_frame_index, s_force& force);

    void advection(double* Ux[], double* Uy[], int frame_index);
    void viscous_diffusion(double* Ux[], double* Uy[], double* scratch, int frame_index);
    void apply_force(double* Ux[], double* Uy[], int frame_index, s_force& force);
    void compute_pressure(double* Ux[], double* Uy[], double* p[], double* scratch, int frame_index);
    void subtract_pressure_gradient(double* Ux[], double* Uy[], double* p[], int frame_index, int p_frame_index);

	void jacobi_loop(double* U[], double* b, int frame_index, double alpha, double rbeta);
	void divergence(double* Ux[], double* Uy[], double* div[], int frame_index);
};

#endif