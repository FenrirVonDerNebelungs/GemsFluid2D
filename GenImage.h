#pragma once
#ifndef GENIMAGE_H
#define GENIMAGE_H
#ifndef BASE_H
#include "Base.h"
#endif
#include "framework.h"

struct s_rgb {
	unsigned char c[3];
};
class GenImage
{
public:
	GenImage(int width, int height, int blow_factor=1, int blow_grid_thickness=0);
	~GenImage();

	s_WH getWidthHeight() { s_WH wh; wh.width = m_width; wh.height = m_height; return wh; }
	BITMAPINFO* getBitmapInfo();
	unsigned char* getImageData() { return m_pImageData; }
	bool makeCopyOfImageData(unsigned char* destBuffer, int destBufferPixelSize=-1/*if >= 0 then checks for correct size in pixels*/);/*copies image data to destBuffer*/
	bool copyImageDataFromBuffer(const unsigned char* srcBuffer, int srcBufferPixelSize = -1/*if >= 0 then checks for correct size in pixels*/);/*copies image data from srcBuffer to m_pImageData*/

	void genTestImage();
	void genNormalizedImage(const double* data);
	void genScaledImage(const double* data, double scale);

	void drawGradientLine(int x1, int y1, int x2, int y2, s_rgb color_start, s_rgb color_end);
	void drawSolidBox(int x_center, int y_center, int dim, s_rgb color);
	void drawBlowSingleColorGrid(s_rgb color_foreground, s_rgb color_background);

	int getBlowFactor() { return m_blow_factor; }
	int getBlowGridThickness() { return m_blow_grid_thickness; }

	s_rgb getColor(double value, double supremum);
private:
	BITMAPINFO* m_pbmi;
	int m_width;
	int m_height;
	unsigned char* m_pImageData;
	int m_blow_factor;
	int m_blow_grid_thickness;

	double getDataSupremum(const double* data);

	s_rgb getColorRB(double value, double supremum);

	/*init utility*/
	void initBitMapHeader(BITMAPINFO* pbmi, int width, int height);
	/*utility*/
	bool isInBlowGridLine(int x, int y);/*returns true if pixel is in grid line where grid line is almost centered about its thickness*/
	bool isInImage(int x, int y);
};
#endif
