#pragma once
#ifndef TEST_H
#define TEST_H
#ifndef CUDAWRAP_H
#include "CUDAWrap.h"
#endif

class Test {
public:
	Test(int num_display_frames=10);
	~Test();

	int runCUDA(double* Ux, double* Uy, double* pressure, s_force& force, int sim_frames);
	
private:
	int m_num_display_frames;
	CUDAWrap* m_pCUDA_wrap;
	double* m_Ux;
	double* m_Uy;
	double* m_p;
	double** m_display_frames;
	void runTestFrame(double* Ux[], double* Uy[], double* p[], double* scratch, int& frame_index, int& p_frame_index, s_force& force);
	void runJacobiTest(double* U[], double* b, int frame_index, double alpha, double rbeta);
};
#endif