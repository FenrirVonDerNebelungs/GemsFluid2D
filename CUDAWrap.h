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
        int blocks_side_dim = 16,//8,   //must be a multiple of 2
        int threads_side_dim = 16, //should be a multiple of 32 but memory overloading, must be a multiple of 2
        double in_delta_t = 1e-3,
        double in_delta_x = 1e-3,
        double in_nu = 1.0e-6, /*max viscosity slow down for walls, about 1/10th of step velocity distance */
        int jacobi_minBlocks_side_dim=2, /*must be such that 2^some power * minBlocks_side dim = blocks_side_dim */
		int jacobi_minThreads_side_dim = 4, /*must be such that 2^some power * minThreads_side dim = threads_side_dim */
        int in_max_jacobi_loops = 200,//50,
        int in_max_jacobi_force_loops = 200,//100,
        double dye_Drho = 0.1
    );
    ~CUDAWrap();

    int runCUDA(double* Ux, double* Uy, double* pressure, double* dye, s_force& force, int sim_frames, double jacobi_filter[]);

	s_WH getGridWidthHeight() { s_WH wh; wh.width = grid_width; wh.height = grid_height; return wh; }

	friend class Test;
protected:

    double* m_dev_Ux[2];/* 2 frames for in and out, each frame is grid_width*grid_height doubles */
    double* m_dev_Uy[2];/* 2 frames for in and out, each frame is grid_width*grid_height doubles */
    double* m_dev_p0[2];/* 2 pressure frames will hold the pressure that counters the advection */
    double* m_dev_p[2];/* 2 frames for in and out, each frame is grid_width*grid_height doubles holds the pressure that counters forces*/
    double* m_dev_scratch;
    double* m_dev_dye[2];

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
    double m_dye_Drho;

    int* advection_indexes[4]; /*raw, or full indexes of corners around back displaced point*/
    double* advection_dist[4]; /*distances scaled to one between displaced point and nearest corners*/

    int m_max_jacobi_loops;
	int m_max_jacobi_force_loops;
	int m_num_jacobi_loops;
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

    int runNV(double* Ux, double* Uy, double* pressure, double* dye, s_force& force, int sim_frames, double jacobi_filter[]);

	void runFrame(
        int& frame_index, 
        int& p_frame_index, 
        int& p_advection_frame_index, 
		int& dye_frame_index,
        s_force& force);

    void advection(int frame_index, int dye_frame_index);
    void viscous_diffusion(int frame_index);
    void apply_force(int frame_index, int dye_frame_index, s_force& force);
    void compute_pressure(double* p[], int frame_index, int p_frame_index);
    void subtract_pressure_gradient(double* p[], int frame_index, int p_frame_index);



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

    /*test functions*/
    void advection_backtrace(
        double* relPos_i,
        double* relPos_j,
        /*const*/ double* Ux[],
        /*const*/ double* Uy[],
        int frame_index);/* test function */
    void bilinearAprox_scaledFrame(
        double* Ux_scaled, 
        double* Uy_scaled, 
        double* dye_scaled,
        const double* Ux, 
        const double* Uy, 
        const double* dye,
        int scale_factor = 6);

    /*general vector math*/
	void gradient(double* dp_dx, double* dp_dy, double* p[], int p_frame_index);
    void divergence(double* div, const double* Ux, const double* Uy);
    void derivative(double* dX, double* dY, const double* Ux, const double* Uy);
    void laplacian(double* lap, double* p[], double* scratch_x, double* scratch_y, int p_frame_index);

    /* init/release  */
    bool initDevMem(int size, double jacobi_filter[]);
	bool copyHostToDeviceMem(int size, double* Ux, double* Uy, double* pressure, double* dye);
	bool copyDeviceToHostMem(
        int size, 
        double* Ux, 
        double* Uy, 
        double* pressure, 
        double* dye, 
        int frame_index, 
        int p_frame_index,
        int dye_frame_index);
	void releaseDevMem();

    /*utility*/
    void add_full_grid(double* dest, const double* src);
    void copyMemory_for_standard_grid(double* dest, const double* src);
    int findLog_base2(int val);
    int find2Pow(int pow);
	int find_stack_WH_and_redFactor(s_WH& wh, int stack_level, int stack_height);
	void find_stack_BlocksNThreads(int& numBlocks_s, int& numThreads_s, int reduction_factor);
};

#endif