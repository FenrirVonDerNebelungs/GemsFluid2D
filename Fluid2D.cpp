#include "Fluid2D.h"

Fluid2D::Fluid2D(
	double mouse_max_delta,
	double mouse_min_delta,
	int sim_frames,
	int num_sims,
	int max_force_frame_duration,
	int force_decay_frames,
	int force_decay_factor,
	double max_allowed_force,
	int max_dye_frames_duration,
	double in_delta_t,
	double in_delta_x,
	double in_nu,
	int blocks_side_dim,
	int threads_side_dim,
	int jacobi_minBlocks_side_dim,
	int jacobi_minThreads_side_dim,
	int in_max_jacobi_loops,
	int in_max_jacobi_force_loops,
	int filter_sigma,
	double dye_intensity
) : 
	m_CUDA_wrap(
		blocks_side_dim, 
		threads_side_dim, 
		in_delta_t, 
		in_delta_x, 
		in_nu, 
		jacobi_minBlocks_side_dim, 
		jacobi_minThreads_side_dim, 
		in_max_jacobi_loops, 
		in_max_jacobi_force_loops),
	m_fluid_animate(max_force_frame_duration, force_decay_frames, force_decay_factor, max_allowed_force, max_dye_frames_duration),
	m_filter(filter_sigma),
	m_file_output_initalized(false),
	m_sim_cnt(0),
	m_mouse_max_delta(mouse_max_delta),
	m_mouse_min_delta(mouse_min_delta),
	m_dye_intensity(dye_intensity),
	m_draw_dye(false),
	m_add_dye_with_force(true),
	m_add_dye_without_force(false),
	m_sim_frames(sim_frames),
	m_num_sims(num_sims),
	m_mouse_clicks(0)
{
	m_images_p = nullptr;
	m_images_dye = nullptr;
	m_Ux = nullptr;
	m_Uy = nullptr;
	m_p = nullptr;
	m_dye = nullptr;
	/*check values*/
	if (sim_frames < 2)
		m_sim_frames = 9;
	if (num_sims < 2)
		m_num_sims = 12;
	/**************/
	s_WH grid_wh = m_CUDA_wrap.getGridWidthHeight();
	m_grid_width = grid_wh.width;
	m_grid_height = grid_wh.height;
	m_size = m_grid_height * m_grid_width;
	if (m_num_sims > 1 && m_grid_width>2 && m_grid_height>2) {
		m_images_p = new GenImage * [m_num_sims];
		m_images_dye = new GenImage * [m_num_sims];
		for (int i = 0; i < m_num_sims; i++) {
			m_images_p[i] = new GenImage(m_grid_width,m_grid_height);
			m_images_dye[i] = new GenImage(m_grid_width,m_grid_height,Gem);
		}
	}
	m_Ux = new double[m_size];
	m_Uy = new double[m_size];
	m_p = new double[m_size];
	m_dye = new double[m_size];
	m_fluid_animate.init(grid_wh);
	std::memset(m_Ux, 0, m_size * sizeof(double));
	std::memset(m_Uy, 0, m_size * sizeof(double));
	std::memset(m_p, 0, m_size * sizeof(double));
	std::memset(m_dye, 0, m_size * sizeof(double));
}
Fluid2D::~Fluid2D() {
	if(m_dye!=nullptr)
		delete[] m_dye;
	if (m_p != nullptr)
		delete[]m_p;
	if (m_Uy != nullptr)
		delete[]m_Uy;
	if (m_Ux != nullptr)
		delete[] m_Ux;
	if(m_images_dye!=nullptr){
		for (int i = 0; i < m_num_sims; i++) {
			if(m_images_dye[i]!=nullptr)
				delete m_images_dye[i];
		}
		delete[] m_images_dye;
	}
	if(m_images_p!=nullptr){
		for (int i = 0; i < m_num_sims; i++) {
			if(m_images_p[i]!=nullptr)
				delete m_images_p[i];
		}
		delete[] m_images_p;
	}
}


bool Fluid2D::initFileOutput(const char* filename) {
	if(m_file_output_initalized)
		return false;
	m_file_output_initalized = true;
	bool initOK = false;
	int num_cache_headers = 4 * m_num_sims;
	int cache_data_len = num_cache_headers * m_size;
	if(filename==nullptr)
		initOK = m_pyTrans.init("Dat/frames.dat", cache_data_len, num_cache_headers);
	else
		initOK = m_pyTrans.init(filename, cache_data_len, num_cache_headers);
	return initOK;
}
void Fluid2D::releaseFileOutput() {
	if(m_file_output_initalized)
		m_pyTrans.release();
	m_file_output_initalized = false;
}
bool Fluid2D::handleMouseSweep(int mouse_end_x, int mouse_end_y_raw, int mouse_start_x, int mouse_start_y_raw) {
	m_mouse_clicks++;
	int mouse_end_y = m_grid_height - mouse_end_y_raw;
	int mouse_start_y = m_grid_height - mouse_start_y_raw;
	if (m_sim_cnt >= m_num_sims)
		m_sim_cnt = 0;/*restarts the simulation for every click*/
	if (mouse_start_x >= (m_grid_width - m_mouse_min_delta) || mouse_start_y >= (m_grid_height - m_mouse_min_delta))
		return false;
	if (mouse_start_x < m_mouse_min_delta || mouse_start_y < m_mouse_min_delta)
		return false;
	int mouse_center_offset_x = mouse_start_x - (m_grid_width / 2);
	int mouse_center_offset_y = mouse_start_y - (m_grid_height / 2);
	double mouse_dx = static_cast<double>(mouse_end_x - mouse_start_x);
	double mouse_dy = static_cast<double>(mouse_end_y - mouse_start_y);
	double dye_intensity = m_add_dye_with_force ? m_dye_intensity : 0.0;
	return applyForce(mouse_dx, mouse_dy, mouse_center_offset_x, mouse_center_offset_y,dye_intensity);
}

bool Fluid2D::runSim() {
	int launchOK = 0;
	double jacobi_filter[g_jacobi_filter_size];
	m_filter.genFilter(jacobi_filter, BASE_JACOBI_EXPANSION_FILTER_HALF_WH);
	if (m_sim_cnt < m_num_sims) {
		std::memset(m_p, 0, m_size * sizeof(double));
		launchOK = m_CUDA_wrap.runCUDA(m_Ux, m_Uy, m_p, m_dye, m_sim_frames, jacobi_filter, &m_fluid_animate);
		if (launchOK == 0) {
			if (m_images_p != nullptr) {
				m_images_p[m_sim_cnt]->genNormalizedImage(m_p);
				m_images_p[m_sim_cnt]->drawMarkerBox(m_sim_cnt, m_num_sims);
			}
			if (m_images_dye != nullptr) {
				m_images_dye[m_sim_cnt]->genNormalizedImage(m_dye);
				m_images_dye[m_sim_cnt]->drawMarkerBox(m_sim_cnt, m_num_sims);
			}
		}
		m_sim_cnt++;
		writeSim();
	}
	else
		return false;
	return launchOK == 0;
}
void Fluid2D::resetSim() {
	std::memset(m_Ux, 0, m_size * sizeof(double));
	std::memset(m_Uy, 0, m_size * sizeof(double));
	std::memset(m_p, 0, m_size * sizeof(double));
	std::memset(m_dye, 0, m_size * sizeof(double));
	m_fluid_animate.setForceActive(false);
	m_fluid_animate.setDyeActive(false);
}
int Fluid2D::writeSim() {
	if(!m_file_output_initalized)
		return 1;
	m_pyTrans.cacheGrid(m_Ux, m_grid_width, m_grid_height, m_sim_cnt, n_PyTrans::U_code, n_PyTrans::X_code, n_PyTrans::end_frame_code);
	m_pyTrans.cacheGrid(m_Uy, m_grid_width, m_grid_height, m_sim_cnt, n_PyTrans::U_code, n_PyTrans::Y_code, n_PyTrans::end_frame_code);
	m_pyTrans.cacheGrid(m_p, m_grid_width, m_grid_height, m_sim_cnt, n_PyTrans::P_code, n_PyTrans::Scalar_code, n_PyTrans::end_frame_code);
	m_pyTrans.cacheGrid(m_dye, m_grid_width, m_grid_height, m_sim_cnt, n_PyTrans::Dye_code, n_PyTrans::Scalar_code, n_PyTrans::end_frame_code);
	m_pyTrans.resetAndWrite();
	return 0;
}
bool Fluid2D::applyForce(double mouse_dx, double mouse_dy, int mouse_center_offset_x, int mouse_center_offset_y, double dye_intensity) {
	double mouse_delta = std::sqrt(mouse_dx * mouse_dx + mouse_dy * mouse_dy);
	if (mouse_delta < m_mouse_min_delta)
		return false;

	if (!m_add_dye_without_force) {
		double mouse_cos = mouse_dx / mouse_delta;
		bool mouse_sin_positive = (mouse_dy >= 0);
		double Force_mag = m_fluid_animate.getMaxAllowedForce();
		if (mouse_delta < m_mouse_max_delta)
			Force_mag *= mouse_delta / m_mouse_max_delta;
		m_fluid_animate.getForce(0, 0, mouse_cos, mouse_sin_positive, Force_mag, mouse_center_offset_x, mouse_center_offset_y);
	}
	else
		m_fluid_animate.getForce(0, 0, 1.0, true, 0.0, mouse_center_offset_x, mouse_center_offset_y);
	bool dye_ok = m_fluid_animate.setDye(0, 0, dye_intensity, mouse_center_offset_x, mouse_center_offset_y);
	return true & dye_ok;
}