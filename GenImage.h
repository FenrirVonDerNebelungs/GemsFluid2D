#pragma once
#ifndef GENIMAGE_H
#define GENIMAGE_H
#ifndef BASE_H
#include "Base.h"
#endif
#include "framework.h"
#include <cmath>

struct s_rgb {
	unsigned char c[3];
};
class GenImage
{
public:
	GenImage(int width, int height);
	~GenImage();

	s_WH getWidthHeight() { s_WH wh; wh.width = m_width; wh.height = m_height; return wh; }
	BITMAPINFO* getBitmapInfo();
	unsigned char* getImageData() { return m_pImageData; }

	void genTestImage();
	void genNormalizedImage(const double* data);
	void genScaledImage(const double* data, double scale);
private:
	BITMAPINFO* m_pbmi;
	int m_width;;
	int m_height;
	unsigned char* m_pImageData;

	double getDataSupremum(const double* data);
	s_rgb getColor(double value, double supremum);
	s_rgb getColorRB(double value, double supremum);
};
#endif
