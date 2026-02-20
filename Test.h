#pragma once
#ifndef TEST_H
#define TEST_H
#ifndef FLUIDANIMATE_H
#include "FluidAnimate.h"
#endif
#ifndef CUDAWRAP_H
#include "CUDAWrap.h"
#endif

class Test {
public:
	Test(int num_display_frames=10);
	~Test();

	int runTest(int sim_frames=1);/*runs the test, returns 0 if successful, otherwise error code*/
	int runCUDA(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames);
	double findErrors();/*assumes runCUDA has already been run, returns max_error*/

	double* getUx() { return m_Ux; }
	double* getUy() { return m_Uy; }
	double* getP() { return m_p; }
	double* getScratch() { return m_scratch; }
	double getCurrentMax() { return m_current_max; }
	double getMaxError() { return m_max_error; }
	double* getDisplayFrame(int frame_index) { return m_display_frames[frame_index]; }
	double* getFrameToDisplay();
	bool preDisplayFrames() { return m_mouse_clicks <= 3; }
	bool switchToDisplayingFrames() { return m_mouse_clicks == 4; }
	bool displayingErrors() { return m_mouse_clicks >= 3+m_num_display_frames && m_mouse_clicks < 3+2*m_num_display_frames; }
	bool switchToDisplayingErrors() { return m_mouse_clicks == 3 + m_num_display_frames; }

	s_WH getGridWidthHeight() { s_WH wh; wh.width = m_pCUDA_wrap->grid_width; wh.height = m_pCUDA_wrap->grid_height; return wh; }
private:
	double m_Ux_max;
	double m_Uy_max;
	double m_p_max;
	double m_current_max;
	int m_mouse_clicks;
	int m_num_display_frames;
	CUDAWrap* m_pCUDA_wrap;
	double* m_Ux;
	double* m_Uy;
	double* m_p;
	double* m_scratch;
	double** m_display_frames;
	double m_max_error;
	FluidAnimate* m_pFluid_animate;
	s_force m_force;

	void runTestFrame(double* Ux[], double* Uy[], double* p[], double* scratch, int& frame_index, int& p_frame_index, s_force& force);
	void runJacobiTest(double* U[], double* b, int frame_index, double alpha, double rbeta);

	void find_max(double& max, const double* data);
	void find_max();
	void set_display_frames_to_errors();/*each display frame is the difference between it and the next, except for the last*/
	void fill_display_frame(int frame_index, const double* dev_data);
	void fill_display_frame(double* pFrame, const double* dev_data);
};
#endif