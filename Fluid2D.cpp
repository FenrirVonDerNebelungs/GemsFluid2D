#include "Fluid2D.h"

Fluid2D::Fluid2D(
	int blocks_side_dim,
	int threads_side_dim,
	int sim_frames
) : 
	m_cuda_wrap(blocks_side_dim, threads_side_dim), 
	m_frame_cnt(0),
	m_blocks_side_dim(blocks_side_dim), 
	m_threads_side_dim(threads_side_dim), 
	m_sim_frames(sim_frames) 
{
	m_Ux = nullptr;
	m_Uy = nullptr;
	m_p = nullptr;
	m_grid_width = blocks_side_dim * threads_side_dim;
	m_grid_height = m_grid_width;
	m_size = m_grid_height * m_grid_width;
	m_Ux = new double[m_size];
	m_Uy = new double[m_size];
	m_p = new double[m_size];
	s_WH grid_wh = { m_grid_width, m_grid_height };
	m_fluid_animate.init(grid_wh);
	std::memset(m_Ux, 0, m_size * sizeof(double));
	std::memset(m_Uy, 0, m_size * sizeof(double));
	std::memset(m_p, 0, m_size * sizeof(double));
}
Fluid2D::~Fluid2D() {
	if (m_p != nullptr)
		delete[]m_p;
	if (m_Uy != nullptr)
		delete[]m_Uy;
	if (m_Ux != nullptr)
		delete[] m_Ux;
}

int Fluid2D::launchCUDA() {
	s_force force = m_fluid_animate.getForce(m_frame_cnt);
	int launchOK = m_cuda_wrap.runCUDA(m_Ux, m_Uy, m_p, force, m_sim_frames);
	if(launchOK==0)
		m_frame_cnt += m_sim_frames;
	return launchOK;
}