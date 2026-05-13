#include "Filter.h"

Filter::Filter(double gaussian_sigma) : m_sigma(gaussian_sigma), m_half_width(0) {
	;
}
Filter::~Filter() {
	;
}

bool Filter::genFilter(double filter[], int target_helf_width) {
	if (target_helf_width < 1 || m_sigma<=0.0)
		return false;
	m_half_width = target_helf_width;
	int width = 2 * target_helf_width;
	double total_weight = 0.0;
	for (int j = 0; j < width; j++) {
		int filter_index = j * width;
		for (int i = 0; i < width; i++) {
			double r = r_from_ij(i, j);
			double cell_weight = getGaussianValue(r);
			filter[filter_index] = cell_weight;
			filter_index++;
			total_weight += cell_weight;
		}
	}
	if (total_weight <= 0.0)
		return false;/*should never happen*/
	int filter_size = width * width;
	for (int i = 0; i < filter_size; i++)
		filter[i] /= total_weight;
	return true;
}
double Filter::r_from_ij(int i, int j) {
	const double max_lo_ij = static_cast<double>(m_half_width - 1);
	double dist_i = static_cast<double>(i) - max_lo_ij;
	dist_i -= 0.5;
	double dist_j = static_cast<double>(j) - max_lo_ij;
	dist_j -= 0.5;
	return sqrt(dist_i * dist_i + dist_j * dist_j);
}
double Filter::getGaussianValue(double r) {
	static const double norm_exp_factor = -1.0 / (2.0 * m_sigma * m_sigma);
	double exp_factor = r * r * norm_exp_factor;
	double gauss_value = exp(exp_factor);
	return gauss_value;
}