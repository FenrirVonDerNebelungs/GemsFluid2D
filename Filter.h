#pragma once
#ifndef FILTER_H
#define FILTER_H

#ifndef BASE_H
#include "Base.h"
#endif

class Filter {
public:
	Filter(double gaussian_sigma=1.0); /*1 is the width of the small cells*/
	~Filter();

	/*generates square filter with gaussian-based weights*/
	bool genFilter(double filter[], int target_half_width);
private:
	double m_sigma;
	double m_half_width;

	double r_from_ij(int i, int j);
	double getGaussianValue(double r);
};

#endif
