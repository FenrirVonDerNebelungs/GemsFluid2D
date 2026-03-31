#pragma once
#ifndef TEST_H
#define TEST_H
#ifndef FLUIDANIMATE_H
#include "FluidAnimate.h"
#endif
#ifndef DRAWVELOCITY_H
#include "drawVelocity.h"
#endif
#ifndef DRAWTEST_H
#include "drawTest.h"
#endif
#ifndef CUDAWRAP_H
#include "CUDAWrap.h"
#endif
#ifndef PYTRANS_H
#include "PyTrans.h"
#endif

class Test {
public:
	Test(int num_display_frames=10, int blow_factor=24);
	~Test();

	int runTest(int sim_frames=10);/*runs the test, returns 0 if successful, otherwise error code*/
	int runCUDA(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames);

	double* getUx() { return m_Ux; }
	double* getUy() { return m_Uy; }
	double* getP() { return m_p; }
	double* getScratch() { return m_scratch; }
	double getCurrentMax() { return m_current_max; }
	double getUnewMax() const { return m_U_new_max; }
	double getMaxError() { return m_max_error; }
	double* getDisplayFrame(int frame_index) { return m_display_frames[frame_index]; }

	GenImage* getTestImage() { return m_drawTest->m_pGenImage; }
	GenImage* handleMouse();
	const WCHAR* getMessage() { return m_drawTest->GetMessage(); }

	s_WH getGridWidthHeight() { s_WH wh; wh.width = m_pCUDA_wrap->grid_width; wh.height = m_pCUDA_wrap->grid_height; return wh; }
	s_WH getBlownWidthHeight() { s_WH wh; wh.width = m_pCUDA_wrap->grid_width * m_blow_factor; wh.height = m_pCUDA_wrap->grid_height * m_blow_factor; return wh; }
private:
	int m_current_frame;
	double m_Ux_max;
	double m_Uy_max;
	double m_U_max;
	double m_div_max;
	double m_Ux_new_max;
	double m_Uy_new_max;
	double m_U_new_max;
	double m_p_max;
	double m_current_max;
	int m_mouse_clicks;
	int m_num_display_frames;
	CUDAWrap* m_pCUDA_wrap;
	int m_blow_factor;
	double* m_Ux;
	double* m_Uy;
	double* m_relPos_x;
	double* m_relPos_y;
	double* m_p;
	double* m_scratch;
	double** m_display_frames;
	double* m_Ux_new;
	double* m_Uy_new;
	double* m_U_div;
	double m_max_error;
	double* m_Ux_bilinear;
	double* m_Uy_bilinear;
	PyTrans* m_pPyTrans;
	FluidAnimate* m_pFluid_animate;
	s_force m_force;
	drawVelocity* m_drawVelocity;
	drawTest* m_drawTest;
	GenImage* m_pCurrent_GenImage;
	double m_image_sup; 

	void runTestFrame(
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
		s_force& force);
	void runTestToPressure(double* Ux[], double* Uy[], double* p[], double* scratch, int& frame_index, int& p_frame_index, s_force& force);
	void runTestAdvection(
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
		s_force& force);
	void runJacobiTest(double* U[], double* b, int frame_index, double alpha, double rbeta, const double* Wx=nullptr, const double* Wy=nullptr);

	void find_max(double& max, const double* data);
	void find_max();
	void set_display_frames_to_errors();/*each display frame is the difference between it and the next, except for the last*/
	void fill_display_frame(int frame_index, const double* dev_data);
	void fill_display_frame(double* pFrame, const double* dev_data);
	void fill_display_frame(double* pFrame, const double* dev_data, int dat_len);

	double* getFrameToDisplay();
	void drawDisplayFrames();
};
#endif