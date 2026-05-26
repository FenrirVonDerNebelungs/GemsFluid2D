#pragma once
#ifndef GENIMAGE_H
#define GENIMAGE_H
#ifndef BASE_H
#include "Base.h"
#endif
#include "framework.h"

enum ColorBarType {
	Rainbow,
	SepiaGlow
};
struct s_rgb {
	unsigned char c[3];
};

class ColorBar {
public:
	ColorBar(ColorBarType type=Rainbow, bool is_centered=false);
	~ColorBar();

	s_rgb getColor(double t_raw);
private:
	s_rgb m_color_marks[6];
	double m_mark_positions[6];/*between 0 and 1 first color must be at 0 last color must be at 1*/
	int m_num_color_marks;
	int m_center_mark_index;
	bool m_is_centered;

	s_rgb m_color_bar[256]; /*rendered color bar*/

	s_rgb getInterpolatedColor(double t);

	s_rgb interpolateColor(s_rgb color1, s_rgb color2, double t);/*t is between 0 and 1*/
	void renderColorBar();/*renders the color bar to m_color_bar for fast lookup*/

	void setRainbowColorBar();
	void setSepiaGlowColorBar();
	void computeHSVfromRGB(const s_rgb& rgb, double& h, double& s, double& v);
	void computeRGBfromHSV(double h, double s, double v, s_rgb& rgb);
};

class GenImage
{
public:
	GenImage(int width, int height, ColorBarType colorBarType=Rainbow, int blow_factor=1, int blow_grid_thickness=0, int marker_box_size=5);
	~GenImage();

	s_WH getWidthHeight() { s_WH wh; wh.width = m_width; wh.height = m_height; return wh; }
	BITMAPINFO* getBitmapInfo();
	unsigned char* getImageData() { return m_pImageData; }
	bool makeCopyOfImageData(unsigned char* destBuffer, int destBufferPixelSize=-1/*if >= 0 then checks for correct size in pixels*/);/*copies image data to destBuffer*/
	bool copyImageDataFromBuffer(const unsigned char* srcBuffer, int srcBufferPixelSize = -1/*if >= 0 then checks for correct size in pixels*/);/*copies image data from srcBuffer to m_pImageData*/

	void genTestImage();
	void genNormalizedImage(const double* data);
	void genScaledImage(const double* data, double min, double max);

	void drawGradientLine(int x1, int y1, int x2, int y2, s_rgb color_start, s_rgb color_end);
	void drawSolidBox(int x_center, int y_center, int dim, s_rgb color);
	void drawMarkerBox(int cnt_current, int cnt_max);

	int getBlowFactor() { return m_blow_factor; }

	s_rgb getColor(double value, double supremum);
	s_rgb getColor(double value, double min, double max);
private:
	BITMAPINFO* m_pbmi;
	ColorBar m_color_bar;
	int m_width;
	int m_height;
	unsigned char* m_pImageData;
	int m_blow_factor;
	int m_marker_box_size;

	double getDataSupremum(const double* data);
	void getDataMinMax(const double* data, double& min, double& max);

	/*init utility*/
	void initBitMapHeader(BITMAPINFO* pbmi, int width, int height);
	/*utility*/
	bool isInImage(int x, int y);

};
#endif
