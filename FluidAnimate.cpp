#include "FluidAnimate.h"
FluidAnimate::FluidAnimate(int max_frame_duration, int decay_frames, double decay_factor, double max_allowed_force, int max_dye_frames_duration) :
	m_grid_width(0),
	m_grid_height(0),
	m_max_frame_duration(max_frame_duration),
	m_decay_frames(decay_frames),
	m_decay_factor(decay_factor),
	m_force_frame_start(0),
	m_force_angle(0.0),
	m_force_cos(0.0),
	m_force_sin(0.0),
	m_force_magnitude(0.0),
	m_max_allowed_force(max_allowed_force),
	m_max_dye_frames_duration(max_dye_frames_duration),
	m_dye_frame_start(0),
	m_dye_intensity(0.0)
{
	m_vars.i = 0;
	m_vars.j = 0;
	m_vars.Fx_c = 0.0;
	m_vars.Fy_c = 0.0;
	m_vars.R = 0.0;
	m_vars.inv_Rsqrd = 0.0;
	m_vars.active = false;
	m_dye.i = 0;
	m_dye.j = 0;
	m_dye.Fx_c = 0.0;
	m_dye.Fy_c = 0.0;
	m_dye.R = 0.0;
	m_dye.inv_Rsqrd = 0.0;
	m_dye.active = false;
}
FluidAnimate::~FluidAnimate() {
	;
}
void FluidAnimate::init(s_WH& grid_wh) {
	m_grid_width = grid_wh.width;
	m_grid_height = grid_wh.height;
}
s_force FluidAnimate::getForce(int frame_cnt, int start_frame, double Ang_cos, bool sin_positive, double Fmag, int center_i, int center_j, double R) {
	m_force_frame_start = start_frame;
	m_force_angle = acos(Ang_cos);
	if (!sin_positive)
		m_force_angle = -m_force_angle;
	m_force_cos = Ang_cos;
	m_force_sin = sqrt(1.0 - Ang_cos * Ang_cos);
	if (!sin_positive)
		m_force_sin = -m_force_sin;
	m_force_magnitude = Fmag;
	m_vars.R = R;
	if(m_vars.R>0.0)
		m_vars.inv_Rsqrd = 1.0 / (R * R);
	else
		m_vars.inv_Rsqrd = 0.0;
	computeForceVars(frame_cnt, center_i, center_j);
	return m_vars;
}
s_force FluidAnimate::updateForce(int frame_cnt) {
	if (m_vars.active && (frame_cnt > (m_max_frame_duration + m_force_frame_start)) ) {
		if(frame_cnt > (m_max_frame_duration + m_decay_frames + m_force_frame_start))
		    m_vars.active = false;
		else {
			m_vars.Fx_c *= m_decay_factor;
			m_vars.Fy_c *= m_decay_factor;
		}
	}
	return m_vars;
}
bool FluidAnimate::setDye(int frame_cnt, int start_frame, double dye_intensity, int center_i, int center_j, double R) {
	m_dye_frame_start = start_frame;
	m_dye_intensity = dye_intensity;
	m_dye.R = R;
	if(m_dye.R>0.0)
		m_dye.inv_Rsqrd = 1.0 / (R * R);
	else
		m_dye.inv_Rsqrd = 0.0;
	return computeDyeVars(frame_cnt, center_i, center_j);
}
s_force FluidAnimate::updateDye(int frame_cnt) {
	if (m_dye.active && (frame_cnt > (m_max_dye_frames_duration + m_dye_frame_start))) {
		m_dye.active = false;
		m_dye.Fx_c = 0.0;
		m_dye.Fy_c = 0.0;
	}
	else if (m_dye.active) {
		m_dye.Fx_c += m_dye_intensity;
		m_dye.Fy_c += m_dye_intensity;
	}
	return m_dye;
}
bool FluidAnimate::computeForceVars(int frame_cnt, int center_offset_i, int center_offset_j) {
	static double grid_width = static_cast<double>(m_grid_width);
	static double grid_height = static_cast<double>(m_grid_height);
	int center_target_i = static_cast<int>(round(grid_width / 2)) + center_offset_i;
	int center_target_j = static_cast<int>(round(grid_height / 2)) + center_offset_j;
	bool center_in_range = inRange_ij_CPU(center_target_i, center_target_j, m_grid_width, m_grid_height);
	int frame_diff = frame_cnt - m_force_frame_start;
	if (center_in_range && frame_diff <= m_max_frame_duration) {
		m_vars.Fx_c = m_force_magnitude * m_force_cos;
		m_vars.Fy_c = m_force_magnitude * m_force_sin;
		m_vars.i = center_target_i;
		m_vars.j = center_target_j;
		m_vars.active = true;
	}
	else {
		m_vars.active = false;
		return false;
	}
	return true;
}
bool FluidAnimate::computeDyeVars(int frame_cnt, int center_offset_i, int center_offset_j) {
	static double grid_width = static_cast<double>(m_grid_width);
	static double grid_height = static_cast<double>(m_grid_height);
	int center_target_i = static_cast<int>(round(grid_width / 2)) + center_offset_i;
	int center_target_j = static_cast<int>(round(grid_height / 2)) + center_offset_j;
	bool center_in_range = inRange_ij_CPU(center_target_i, center_target_j, m_grid_width, m_grid_height);
	int frame_diff = frame_cnt - m_force_frame_start;
	if (center_in_range && frame_diff <= m_max_frame_duration) {
		m_dye.Fx_c = m_dye_intensity;
		m_dye.Fy_c = m_dye_intensity;
		m_dye.i = center_target_i;
		m_dye.j = center_target_j;
		m_dye.active = true;
	}
	else {
		m_dye.active = false;
		return false;
	}
	return true;
}

