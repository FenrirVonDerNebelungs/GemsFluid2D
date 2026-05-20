#pragma once
#ifndef TEST_H
#define TEST_H
#ifndef FLUIDANIMATE_H
#include "../FluidAnimate.h"
#endif
#ifndef DRAWVELOCITY_H
#include "drawVelocity.h"
#endif
#ifndef DRAWTEST_H
#include "drawTest.h"
#endif
#ifndef FILTER_H
#include "../Filter.h"
#endif
#ifndef CUDAWRAP_H
#include "../CUDAWrap.h"
#endif
#ifndef PYTRANS_H
#include "../PyTrans.h"
#endif

class Test {
public:
	Test(int num_display_frames=10, int blow_factor=2, int loop_modulus_divider=50);
	~Test();

	int runTest(int sim_frames=10);/*runs the test, returns 0 if successful, otherwise error code*/
	int runCUDA(double* Ux, double* Uy, double* pressure, double* pressure0, int sim_frames, double jacobi_filter[]);

	double* getUx() { return m_Ux; }
	double* getUy() { return m_Uy; }
	double* getP() { return m_p; }

	GenImage* handleMouse();

	s_WH getGridWidthHeight() { s_WH wh; wh.width = m_pCUDA_wrap->grid_width; wh.height = m_pCUDA_wrap->grid_height; return wh; }
	s_WH getBlownWidthHeight() { s_WH wh; wh.width = m_pCUDA_wrap->grid_width * m_blow_factor; wh.height = m_pCUDA_wrap->grid_height * m_blow_factor; return wh; }

	GenImage* getTestImage() { return m_pDraw->getGenImage(); }
	const WCHAR* getMessage() { return m_pDraw->GetMessage(); }
private:
	int m_current_frame;
	double m_Ux_max;
	double m_Uy_max;
	double m_U_max;
	double m_p_max;
	int m_mouse_clicks;
	CUDAWrap* m_pCUDA_wrap;
	int m_blow_factor;
	double* m_Ux;
	double* m_Uy;
	double* m_p;
	double* m_dye;
	double* m_Ux_bilinear;
	double* m_Uy_bilinear;
	double* m_dye_bilinear;
	PyTrans* m_pPyTrans;
	FluidAnimate* m_pFluid_animate;
	Filter m_filter;
	s_force m_force;

	drawTest* m_pDraw;
	GenImage* m_pCurrent_GenImage;
	double    m_current_max;

	double* m_host_scratch[4];
	int* m_host_scratch_int;
	int m_loop_modulus_divider;

	void runFrame(
		int& frame_index, 
		int& p_frame_index,
		int& p_advection_frame_index,
		int& dye_frame_index,
		s_force& force,
		s_force& dye_brush);
	void runTestAdvectionFrame(
		int& frame_index,
		int& p_frame_index,
		int& p_advection_frame_index,
		int& dye_frame_index,
		s_force& force,
		s_force& dye_brush);
	void compute_pressure(double* p[], int frame_index, int p_frame_index);
	void jacobi_loop(
		double* X[], /*frame in is X[0], frame out is also X[0], X[1] in is swap, the pointer X[0] may = the oringinal pointer X[1]*/
		const double* b,
		double alpha,
		double rbeta,
		int numBlocks_s,
		int numThreads_s,
		const double* Wx = nullptr,
		const double* Wy = nullptr,
		int stack_index=0);
	void jacobi_run(
		double* X[],
		const double* b,
		const double* Wx,
		const double* Wy,
		int frame_index);
	void send_jacobi_stacks_to_cache(int jacobi_frame);
	int getCacheLen(int& num_headers);

	void find_max(double& max, const double* data);
	void find_max();
	void send_frame_to_host(double* pFrame, const double* dev_data);
	void send_frame_to_host(double* pFrame, const double* dev_data, int dat_len);
	void send_frame_to_host(int* pFrame, const int* dev_data);

	double* getFrameToDisplay();
	void drawDisplayFrames();
};
#endif