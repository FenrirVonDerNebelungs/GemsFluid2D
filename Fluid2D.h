#pragma once
#ifndef FLUID2D_H
#define FLUID2D_H
#ifndef FLUIDANIMATE_H
#include "FluidAnimate.h"
#endif
#ifndef CUDAWRAP_H
#include "CUDAWrap.h"
#endif

class Fluid2D {
public:
	Fluid2D(
		int blocks_side_dim = 8,
		int threads_side_dim = 16,
		int sim_frames = 3
	);
	~Fluid2D();

	int launchCUDA();

	s_WH getGridWidthHeight() { s_WH wh; wh.width = m_grid_width; wh.height = m_grid_height; return wh; }
	double* getUx() { return m_Ux; }
	double* getUy() { return m_Uy; }
	double* getP() { return m_p; }
	int getFrameCount() { return m_frame_cnt; }
private:
	CUDAWrap m_cuda_wrap;
	FluidAnimate m_fluid_animate;
	int m_frame_cnt;
	double* m_Ux;
	double* m_Uy;
	double* m_p;
	int m_blocks_side_dim;
	int m_threads_side_dim;
	int m_sim_frames;
	int m_grid_width;
	int m_grid_height;
	int m_size;

};
#endif
