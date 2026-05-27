#include "GenImage.h"

ColorBar::ColorBar(ColorBarType type, bool is_centered) : m_is_centered(is_centered) {
	switch (type) {
	case Rainbow:
		setRainbowColorBar();
		break;
	case SepiaGlow:
		setSepiaGlowColorBar();
		break;
	case Gem:
		setGemColorBar();
		break;
	default:
		setRainbowColorBar();
		break;
	}
}
ColorBar::~ColorBar() {
	;
}
s_rgb ColorBar::getColor(double t_raw) {
	double t = m_is_centered ? (1.0 + t_raw) / 2.0 : t_raw;
	int index = static_cast<int>(t * 255.0);
	if (index < 0)
		index = 0;
	if (index > 255)
		index = 255;
	return m_color_bar[index];
}
s_rgb ColorBar::getInterpolatedColor(double t) {
	if (t < 0.0)
		t = 0.0;
	if (t <= m_mark_positions[0])
		return m_color_marks[0];
	if (t >= m_mark_positions[m_num_color_marks - 1])
		return m_color_marks[m_num_color_marks - 1];
	for (int i = 0; i < m_num_color_marks - 1; i++) {
		if (t > m_mark_positions[i] && t <= m_mark_positions[i + 1]) {
			double local_t = (t - m_mark_positions[i]) / (m_mark_positions[i + 1] - m_mark_positions[i]);
			return interpolateColor(m_color_marks[i], m_color_marks[i + 1], local_t);
		}
	}
	return { 0, 0, 0 }; // Default to black if something goes wrong
}
s_rgb ColorBar::interpolateColor(s_rgb color1, s_rgb color2, double t) {
	s_rgb result;
	/*rotate color while preserving luminance*/
	double h1, s1, v1, h2, s2, v2;
	computeHSVfromRGB(color1, h1, s1, v1);
	computeHSVfromRGB(color2, h2, s2, v2);
	double h_diff = h2 - h1;
	double new_h = h1 + t * h_diff;
	double new_s = s1 + t * (s2 - s1);
	double new_v = v1 + t * (v2 - v1);
	computeRGBfromHSV(new_h, new_s, new_v, result);
	return result;
}
void ColorBar::renderColorBar() {
	for (int i = 0; i < 256; i++) {
		double t_raw = static_cast<double>(i) / 255.0;
		m_color_bar[i] = getInterpolatedColor(t_raw);
	}
}
void ColorBar::reverseColorBar() {
	for (int i = 0; i < 256; i++) {
		char c0 = m_color_bar[i].c[0];
		char c2 = m_color_bar[i].c[2];
		m_color_bar[i].c[0] = c2;
		m_color_bar[i].c[2] = c0;
	}
}
void ColorBar::setRainbowColorBar() {
	m_color_marks[0] = { 200,0,0 };//{ 255, 0, 0 }; // Red
	m_color_marks[1] = { 255, 100,0 };//{ 255, 127, 0 }; // Orange
	m_color_marks[2] = { 255, 255, 0 }; // Yellow
	m_color_marks[3] = { 0,200,95 };//{ 0, 255, 0 }; // Green
	m_color_marks[4] = { 0, 0, 255 }; // Blue
	m_mark_positions[0] = 0.0;
	m_mark_positions[1] = 0.25;
	m_mark_positions[2] = 0.5;
	m_mark_positions[3] = 0.75;
	m_mark_positions[4] = 1.0;
	m_num_color_marks = 5;
	m_center_mark_index = 2;
	renderColorBar();
}
void ColorBar::setSepiaGlowColorBar() {
	m_color_marks[5] = { 255,255,224 };//{ 255, 255, 224 }; // Light Yellow
	m_color_marks[4] = { 255,228,196 };//{ 255, 228, 196 }; // Bisque
	m_color_marks[3] = { 200, 192, 185 };//{ 255, 218, 185 }; // Peach Puff
	m_color_marks[2] = { 130, 120, 110 };//{ 255, 160, 122 }; // Light Salmon
	m_color_marks[1] = {95, 100, 105 };//{ 255, 127, 80 }; // Coral
	m_color_marks[0] = { 0, 10, 30 };//{ 139, 69, 19 }; // Saddle Brown
	m_mark_positions[0] = 0.0;
	m_mark_positions[1] = 0.2;
	m_mark_positions[2] = 0.4;
	m_mark_positions[3] = 0.6;
	m_mark_positions[4] = 0.8;
	m_mark_positions[5] = 1.0;
	m_num_color_marks = 6;
	m_center_mark_index = 0;
	renderColorBar();
}
void ColorBar::setGemColorBar() {
	m_color_marks[5] = { 215,210,225 };//{ 255, 255, 224 }; // Light Yellow
	m_color_marks[4] = { 155,121,225 };//{ 255, 228, 196 }; // Bisque
	m_color_marks[3] = { 100, 90, 200 };//{ 255, 218, 185 }; // Peach Puff
	m_color_marks[2] = { 100, 46, 137};//{ 255, 160, 122 }; // Light Salmon
	m_color_marks[1] = { 76, 6, 86 };//{ 126, 39, 60 }; // Coral
	m_color_marks[0] = { 0, 0, 0 };//{ 139, 69, 19 }; // Saddle Brown
	m_mark_positions[0] = 0.0;
	m_mark_positions[1] = 0.2;
	m_mark_positions[2] = 0.4;
	m_mark_positions[3] = 0.6;
	m_mark_positions[4] = 0.8;
	m_mark_positions[5] = 1.0;
	m_num_color_marks = 6;
	m_center_mark_index = 0;
	renderColorBar();
	reverseColorBar();
}
void ColorBar::computeHSVfromRGB(const s_rgb& rgb, double& h, double& s, double& v) {
	double r = static_cast<double>(rgb.c[0]) / 255.0;
	double g = static_cast<double>(rgb.c[1]) / 255.0;
	double b = static_cast<double>(rgb.c[2]) / 255.0;
	double max = r>g ? r: g;
	max = max > b ? max : b;
	double min = r < g ? r : g;
	min = min < b ? min : b;
	v = max;
	double delta = max - min;
	s = (max == 0) ? 0 : (delta / max);
	if (delta == 0) {
		h = 0; // Undefined hue
	} else {
		if (max == r) {
			h = 60 * (g - b) / delta + 360;
			h = fmod(h, 360);
		} else if (max == g) {
			h = 60 * (b - r) / delta + 120;
			h = fmod(h, 360);
		} else { // max == b
			h = 60 * (r - g) / delta + 240;
			h = fmod(h, 360);
		}
	}
}
void ColorBar::computeRGBfromHSV(double h, double s, double v, s_rgb& rgb) {
	double c = v * s;
	double x = c * (1 - std::fabs(fmod(h / 60.0, 2) - 1));
	double m = v - c;
	double r_prime, g_prime, b_prime;
	if (h < 60) {
		r_prime = c; g_prime = x; b_prime = 0;
	} else if (h < 120) {
		r_prime = x; g_prime = c; b_prime = 0;
	} else if (h < 180) {
		r_prime = 0; g_prime = c; b_prime = x;
	} else if (h < 240) {
		r_prime = 0; g_prime = x; b_prime = c;
	} else if (h < 300) {
		r_prime = x; g_prime = 0; b_prime = c;
	} else {
		r_prime = c; g_prime = 0; b_prime = x;
	}
	rgb.c[0] = static_cast<unsigned char>(std::round((r_prime + m) * 255));
	rgb.c[1] = static_cast<unsigned char>(std::round((g_prime + m) * 255));
	rgb.c[2] = static_cast<unsigned char>(std::round((b_prime + m) * 255));
}
GenImage::GenImage(int width, int height, ColorBarType colorBarType, int blow_factor, int blow_grid_thickness, int marker_box_size) : 
	m_color_bar(colorBarType),
	m_marker_box_size(marker_box_size)
{
	m_width = width;
	m_height = height;
	int nBitsPerPixel = 24;
	//int nBytesPerPixel = nBitsPerPixel / 8;
	//int nBytesPerLine = ((m_width * nBytesPerPixel + 3) / 4) * 4; 
	//int nImageSize = nBytesPerLine * m_height;
	m_pbmi = reinterpret_cast<BITMAPINFO*>( GlobalAlloc(GPTR, sizeof(BITMAPINFOHEADER)) ); // +sizeof(RGBQUAD) * colorsUsed);  GPTR is GMEM_FIXED | GMEM_ZEROINIT which means bytes are in fixed memory referenced by a pointer and mem is zeroed
	initBitMapHeader(m_pbmi, m_width, m_height);
	m_pImageData = new unsigned char[m_width * m_height * nBitsPerPixel / 8];

	m_blow_factor = blow_factor;
}
void GenImage::initBitMapHeader(BITMAPINFO* pbmi, int width, int height) {
	int nBitsPerPixel = 24;
	//int nBytesPerPixel = nBitsPerPixel / 8;
	//int nBytesPerLine = ((m_width * nBytesPerPixel + 3) / 4) * 4; 
	//int nImageSize = nBytesPerLine * m_height;

	pbmi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	pbmi->bmiHeader.biWidth = width;
	pbmi->bmiHeader.biHeight = height;//-height;  autcomplete made this negative and added the comment "top-down image" see what happens with this positive
	pbmi->bmiHeader.biPlanes = 1;
	pbmi->bmiHeader.biBitCount = nBitsPerPixel;
	pbmi->bmiHeader.biCompression = BI_RGB;/* no compression*/
	pbmi->bmiHeader.biSizeImage = 0;/* can be 0 for BI_RGB (not compressed) bitmaps*/
	pbmi->bmiHeader.biXPelsPerMeter = 0;
	pbmi->bmiHeader.biYPelsPerMeter = 0;
	pbmi->bmiHeader.biClrUsed = 0;/* we are using all 16 million colors; use the maximum number of colors for the bitdepth 2^32 */
	pbmi->bmiHeader.biClrImportant = 0;
}
GenImage::~GenImage()
{
	if(m_pImageData!=nullptr)
	{
		delete[] m_pImageData;
		m_pImageData = nullptr;
	}
	if (m_pbmi != nullptr)
	{
		GlobalFree(m_pbmi);
		m_pbmi = nullptr;
	}

}
BITMAPINFO* GenImage::getBitmapInfo()
{
	return m_pbmi;
}
bool GenImage::makeCopyOfImageData(unsigned char* destBuffer, int destBufferPixelSize)
{
	if(destBuffer == nullptr)
		return false; // Invalid buffer
	if (destBufferPixelSize!=-1 && destBufferPixelSize < m_width * m_height)
		return false;
	size_t requiredSize = static_cast<size_t>(m_width * m_height * 3);// 3 bytes per pixel for 24-bit RGB
	std::memcpy(destBuffer, m_pImageData, requiredSize);
	return true;
}
bool GenImage::copyImageDataFromBuffer(const unsigned char* srcBuffer, int srcBufferPixelSize)
{
	if (srcBuffer == nullptr)
		return false; // Invalid buffer
	if (srcBufferPixelSize!=-1 && srcBufferPixelSize < m_width * m_height)
		return false;
	size_t requiredSize= static_cast<size_t>(m_width * m_height * 3);// 3 bytes per pixel for 24-bit RGB
	std::memcpy(m_pImageData, srcBuffer, requiredSize);
	return true;
}
void GenImage::genTestImage()
{
	for (int y = 0; y < m_height; y++)
	{
		for (int x = 0; x < m_width; x++)
		{
			int index = (y * m_width + x) * 3;
			double x_Norm = static_cast<double>(x) / static_cast<double>(m_width);
			double y_Norm = static_cast<double>(y) / static_cast<double>(m_height);	
			m_pImageData[index] = static_cast<unsigned char>(std::round(255.0 * x_Norm)); // Red
			m_pImageData[index + 1] = 0; // Green
			m_pImageData[index + 2] = static_cast<unsigned char>(std::round(255.0 * y_Norm)); // Blue
		}
	}
}
void GenImage::genNormalizedImage(const double* data) {
	if (data == nullptr)
		return;
	double max, min;
	getDataMinMax(data, min, max);
	genScaledImage(data, min, max);
}
void GenImage::genScaledImage(const double* data, double  min, double max) {
	for (int y = 0; y < m_height; y++) {
		for (int x = 0; x < m_width; x++) {
			int index = (y * m_width + x) * 3;
			double value = data[y * m_width + x];
			s_rgb color = getColor(value, min, max);
			m_pImageData[index] = color.c[0]; // Red, actually Blue
			m_pImageData[index + 1] = color.c[1]; // Green
			m_pImageData[index + 2] = color.c[2]; // Blue, actually Red
		}
	}
}
double GenImage::getDataSupremum(const double* data) {
	double supremum = 0.0;
	int size = m_width * m_height;
	for (int i = 0; i < size; i++) {
		double abs_data = std::abs(data[i]);
		if (abs_data > supremum) {
			supremum = abs_data;
		}
	}
	return supremum;
}
void GenImage::getDataMinMax(const double* data, double& min, double& max) {
	double supremum = getDataSupremum(data);
	min = supremum;
	max = -supremum;
	int size = m_width * m_height;
	for (int i = 0; i < size; i++) {
		double abs_data = std::abs(data[i]);
		if (abs_data > supremum) {
			supremum = abs_data;
		}
		if (data[i] < min) {
			min = data[i];
		}
		if (data[i] > max) {
			max = data[i];
		}
	}
	return;
}
s_rgb GenImage::getColor(double value, double supremum) {
	double color_t = value / supremum;
	return m_color_bar.getColor(color_t);
}
s_rgb GenImage::getColor(double value, double min, double max) {
	double color_t = (value - min) / (max - min);
	return m_color_bar.getColor(color_t);
}

void GenImage::drawGradientLine(int x1, int y1, int x2, int y2, s_rgb color_start, s_rgb color_end) {
	if (!isInImage(x1, y1) || !isInImage(x2, y2))
		return;
	int dx = x2 - x1;
	int dy = y2 - y1;
	unsigned int ab_dx = std::abs(dx);
	unsigned int ab_dy = std::abs(dy);
	int steps = (ab_dx>ab_dy) ? ab_dx : ab_dy;
	if (steps == 0) {
		return; // Avoid division by zero
	}
	double x_inc = static_cast<double>(dx) / steps;
	double y_inc = static_cast<double>(dy) / steps;
	for (int i = 0; i <= steps; i++) {
		int x = static_cast<int>(std::round(x1 + i * x_inc));
		int y = static_cast<int>(std::round(y1 + i * y_inc));
		double t = static_cast<double>(i) / steps;
		s_rgb color;
		color.c[0] = static_cast<unsigned char>(std::round((1 - t) * color_start.c[0] + t * color_end.c[0]));
		color.c[1] = static_cast<unsigned char>(std::round((1 - t) * color_start.c[1] + t * color_end.c[1]));
		color.c[2] = static_cast<unsigned char>(std::round((1 - t) * color_start.c[2] + t * color_end.c[2]));
		int index = (y * m_width + x) * 3;
		m_pImageData[index] = color.c[0]; // Red
		m_pImageData[index + 1] = color.c[1]; // Green
		m_pImageData[index + 2] = color.c[2]; // Blue
	}
}
void GenImage::drawSolidBox(int x_center, int y_center, int dim, s_rgb color) {
	int half_dim = dim / 2;
	for (int y = y_center - half_dim; y <= y_center + half_dim; y++) {
		for (int x = x_center - half_dim; x <= x_center + half_dim; x++) {
			if (x >= 0 && x < m_width && y >= 0 && y < m_height) {
				int index = (y * m_width + x) * 3;
				m_pImageData[index] = color.c[0]; // Red
				m_pImageData[index + 1] = color.c[1]; // Green
				m_pImageData[index + 2] = color.c[2]; // Blue
			}
		}
	}
}
void GenImage::drawMarkerBox(int cnt_current, int cnt_max) {
	double sup = static_cast<double>(cnt_max);
	double val = static_cast<double>(cnt_current);
	s_rgb mark_col = getColor(val, sup);
	int marker_offset = (m_marker_box_size / 2 + m_marker_box_size % 2);
	if (cnt_current == (cnt_max-1)) {
		mark_col.c[0] = 255;
		mark_col.c[1] = 120;
		mark_col.c[2] = 255;
	}
	drawSolidBox(marker_offset, marker_offset, m_marker_box_size, mark_col);
}
bool GenImage::isInImage(int x, int y) {
	if (x < 0 || y < 0)
		return false;
	if (x >= m_width || y >= m_height)
		return false;
	return true;
}