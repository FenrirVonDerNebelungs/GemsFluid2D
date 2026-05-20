#pragma once
#ifndef FLUID2D_H
#define FLUID2D_H
#ifndef FILTER_H
#include "Filter.h"
#endif
#ifndef CUDAWRAP_H
#include "CUDAWrap.h"
#endif
#ifndef GENIMAGE_H
#include "GenImage.h"
#endif
#ifndef PYTRANS_H
#include "PyTrans.h"
#endif

class Fluid2D {
public:
	Fluid2D(
		double mouse_max_delta=20,
		double mouse_min_delta=5,
		int sim_frames = 20,
		int num_sims=5,
		int max_force_frame_duration = 6,
		int force_decay_frames = 3,
		int force_decay_factor = 0.1,
		int max_dye_frames_duration = 3,
		double in_delta_t = 1.0e-3,
		double in_delta_x = 1.0e-3,
		double in_nu = 1.0e-4,
		int blocks_side_dim = 16,
		int threads_side_dim = 16,
		int jacobi_minBlocks_side_dim = 2,
		int jacobi_minThreads_side_dim = 4,
		int in_max_jacobi_loops=50,
		int in_max_jacobi_force_loops=100,
		int filter_sigma = 1.0
	);
	~Fluid2D();

	int launchCUDA();

	bool applyForce(int mouse_end_x=0, int mouse_end_y=0, int mouse_center_x=-1, int mouse_center_y=-1, double dye_intensity=0.1);
	bool initFileOutput(const char* filename=nullptr);
	void releaseFileOutput();
	GenImage* handleMouse();

	s_WH getGridWidthHeight() { s_WH wh; wh.width = m_grid_width; wh.height = m_grid_height; return wh; }
	double* getUx() { return m_Ux; }
	double* getUy() { return m_Uy; }
	double* getP() { return m_p; }
	double* getDye() { return m_dye; }

	GenImage* getImagePtr_p(int sim_frame = 0) { return m_images_p != nullptr ? m_images_p[sim_frame] : nullptr; }
	GenImage* getImagePtr_dye(int sim_frame = 0) { return m_images_dye != nullptr ? m_images_dye[sim_frame] : nullptr; }
private:
	CUDAWrap     m_CUDA_wrap;
	FluidAnimate m_fluid_animate;
	Filter       m_filter;
	GenImage** m_images_p;
	GenImage** m_images_dye;
	PyTrans m_pyTrans;
	int m_sim_cnt;
	double* m_Ux;
	double* m_Uy;
	double* m_p;
	double* m_dye;
	double m_mouse_max_delta;
	double m_mouse_min_delta;
	int m_sim_frames;
	int m_num_sims;
	int m_grid_width;
	int m_grid_height;
	int m_size;

	int m_mouse_clicks;
};
#endif
