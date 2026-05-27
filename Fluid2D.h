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
		double mouse_max_delta=60,
		double mouse_min_delta=5,
		int sim_frames = 9, /*number of loops before image is displayed*/
		int num_sims=12, /*total number of loops in sim is sim_frames * num_sims */
		int max_force_frame_duration = 6,
		int force_decay_frames = 3,
		int force_decay_factor = 0.1,
		double max_allowed_force = 300.0,
		int max_dye_frames_duration = 3,
		double in_delta_t = 1.0e-3,
		double in_delta_x = 1.0e-3,
		double in_nu = 1.0e-5,
		int blocks_side_dim = 16,
		int threads_side_dim = 16,
		int jacobi_minBlocks_side_dim = 2,
		int jacobi_minThreads_side_dim = 4,
		int in_max_jacobi_loops=24,
		int in_max_jacobi_force_loops=64,
		int filter_sigma = 1.0,
		double dye_intensity=0.1
	);
	~Fluid2D();

	bool runSim();/*true if images have been updated*/
	void advanceSim() { if (m_sim_cnt >= m_num_sims) m_sim_cnt = 0; }
	void resetSim();

	bool initFileOutput(const char* filename=nullptr);
	void releaseFileOutput();
	bool handleMouseSweep(int mouse_end_x=0, int mouse_end_y_raw=0, int mouse_start_x=-1, int mouse_start_y_raw=-1);
	void setDisplayDye() { m_draw_dye = true; }
	void setDisplayP() { m_draw_dye = false; }
	void toggleAddDyeWithForce() { m_add_dye_with_force = !m_add_dye_with_force; }
	void toggleAddDyeWithoutForce() { m_add_dye_without_force = !m_add_dye_without_force; }

	GenImage* getCurrentImage() { return m_draw_dye ? getCurrentImage_dye() : getCurrentImage_p(); }
	GenImage* getCurrentImage_p() { return (m_images_p != nullptr && m_sim_cnt>0) ? m_images_p[m_sim_cnt-1] : nullptr; }
	GenImage* getCurrentImage_dye() { return (m_images_dye != nullptr && m_sim_cnt>0) ? m_images_dye[m_sim_cnt-1] : nullptr; }
	s_WH getGridWidthHeight() { s_WH wh; wh.width = m_grid_width; wh.height = m_grid_height; return wh; }
	double* getUx() { return m_Ux; }
	double* getUy() { return m_Uy; }
	double* getP() { return m_p; }
	double* getDye() { return m_dye; }

	void toggleSlidingWallActive() { m_fluid_animate.slidingWallActive() ? m_fluid_animate.setSlidingWallActive(false) : m_fluid_animate.setSlidingWallActive(); }
	bool slidingWallActive() { return m_fluid_animate.slidingWallActive(); }
	s_force getSlidingWallU() { return m_fluid_animate.getSlidingWallU(); }

private:
	CUDAWrap     m_CUDA_wrap;
	FluidAnimate m_fluid_animate;
	Filter       m_filter;
	GenImage** m_images_p;
	GenImage** m_images_dye;
	PyTrans m_pyTrans;
	bool m_file_output_initalized;
	int m_sim_cnt;
	double* m_Ux;
	double* m_Uy;
	double* m_p;
	double* m_dye;
	double m_mouse_max_delta;
	double m_mouse_min_delta;
	double m_dye_intensity;
	bool m_draw_dye;
	bool m_add_dye_with_force;
	bool m_add_dye_without_force;
	int m_sim_frames;
	int m_num_sims;
	int m_grid_width;
	int m_grid_height;
	int m_size;

	int m_mouse_clicks;

	int writeSim();
	bool applyForce(double mouse_dx, double mouse_dy, int mouse_center_offset_x, int mouse_center_offset_y, double dye_intensity = 0.1);
	GenImage* getImagePtr_p(int sim_frame = 0) { return m_images_p != nullptr ? m_images_p[sim_frame] : nullptr; }
	GenImage* getImagePtr_dye(int sim_frame = 0) { return m_images_dye != nullptr ? m_images_dye[sim_frame] : nullptr; }
};
#endif
