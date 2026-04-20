#pragma once
#ifndef FLUIDANIMATE_H
#define FLUIDANIMATE_H
#endif
#ifndef BASE_H
#include "Base.h"
#endif

class FluidAnimate {
public:
	FluidAnimate(int max_frame_duration=3);
	~FluidAnimate();

	void init(s_WH& grid_wh);

	s_force getForce(int frame_cnt, int start_frame=0, double Ang=0.174533, double Fmag=1.0, int center_i=-4, int center_j=3, double R=10.0);//Fmag was originally 1
private:
	int m_grid_width;
	int m_grid_height;
	int m_max_frame_duration;

	int m_force_frame_start;
	double m_force_angle;
	double m_force_magnitude;
	int m_center_offset_i;
	int m_center_offset_j;

	s_force m_vars;

	bool computeForceVars(int frame_cnt);
};