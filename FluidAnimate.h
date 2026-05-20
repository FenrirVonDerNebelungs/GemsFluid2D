#pragma once
#ifndef FLUIDANIMATE_H
#define FLUIDANIMATE_H
#endif
#ifndef BASE_H
#include "Base.h"
#endif

class FluidAnimate {
public:
	FluidAnimate(
		int max_frame_duration=6,
		int decay_frames=3,
		double decay_factor=0.1,
		double max_allowed_force = 300.0,
		int max_dye_frames_duration=3);
	~FluidAnimate();

	void init(s_WH& grid_wh);

	s_force getForce(int frame_cnt, int start_frame=0, double Ang_cos=-0.2, bool sin_positive=false, double Fmag=300.0, int center_i=4, int center_j=10, double R=16.0);//Fmag was originally 1
	s_force updateForce(int frame_cnt);
	double getMaxAllowedForce() { return m_max_allowed_force; }
	double getDyeIntensity() { return m_dye_intensity; }	
	bool setDye(int frame_cnt, int start_frame=0, double dye_intensity=0.1, int center_i=4, int center_j=10, double R=3.0);
	s_force updateDye(int frame_cnt);
private:
	int m_grid_width;
	int m_grid_height;
	int m_max_frame_duration;
	int m_decay_frames;
	double m_decay_factor;

	int m_force_frame_start;
	double m_force_angle;
	double m_force_cos;
	double m_force_sin;
	double m_force_magnitude;
	double m_max_allowed_force;

	s_force m_vars;

	int m_max_dye_frames_duration;
	int m_dye_frame_start;
	double m_dye_intensity;
	s_force m_dye;

	bool computeForceVars(int frame_cnt, int center_offset_i, int center_offset_j);
	bool computeDyeVars(int frame_cnt, int center_offset_i, int center_offset_j);
};