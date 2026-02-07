#include "FluidAnimate.h"
FluidAnimate::FluidAnimate(int max_frame_duration) : m_max_frame_duration(max_frame_duration) {
	m_grid_width = 0;
	m_grid_height = 0;
	m_force_frame_start = 0;
	m_force_angle = 0.0;
	m_force_magnitude = 0.0;
	m_center_offset_i = 0;
	m_center_offset_j = 0;
	m_vars.i = 0;
	m_vars.j = 0;
	m_vars.Fx_c = 0.0;
	m_vars.Fy_c = 0.0;
	m_vars.R = 0.0;
	m_vars.inv_Rsqrd = 0.0;
	m_vars.active = false;
}
FluidAnimate::~FluidAnimate() {
	;
}
void FluidAnimate::init(s_WH& grid_wh) {
	m_grid_width = grid_wh.width;
	m_grid_height = grid_wh.height;
}
s_force FluidAnimate::getForce(int frame_cnt, int start_frame, double Ang, double Fmag, int center_i, int center_j, double R) {
	m_force_frame_start = start_frame;
	m_force_angle = Ang;
	m_force_magnitude = Fmag;
	m_center_offset_i = center_i;
	m_center_offset_j = center_j;
	m_vars.R = R;
	if(m_vars.R>0.0)
		m_vars.inv_Rsqrd = 1.0 / (R * R);
	else
		m_vars.inv_Rsqrd = 0.0;
	computeForceVars(frame_cnt);
	return m_vars;
}
bool FluidAnimate::computeForceVars(int frame_cnt) {
	static double grid_width = static_cast<double>(m_grid_width);
	static double grid_height = static_cast<double>(m_grid_height);
	int center_target_i = static_cast<int>(round(grid_width / 2)) + m_center_offset_i;
	int center_target_j = static_cast<int>(round(grid_height / 2)) + m_center_offset_j;
	bool center_in_range = inRange_ij_CPU(center_target_i, center_target_j, m_grid_width, m_grid_height);
	int frame_diff = frame_cnt - m_force_frame_start;
	if (center_in_range && frame_diff <= m_max_frame_duration) {
		m_vars.Fx_c = m_force_magnitude * std::cos(m_force_angle);
		m_vars.Fy_c = m_force_magnitude * std::sin(m_force_angle);
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