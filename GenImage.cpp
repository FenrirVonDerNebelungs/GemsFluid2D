#include "GenImage.h"

GenImage::GenImage(int width, int height, int blow_factor, int blow_grid_thickness)
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
	m_blow_grid_thickness = blow_grid_thickness;
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
	double supremum = getDataSupremum(data);
	genScaledImage(data, supremum);
}
void GenImage::genScaledImage(const double* data, double scale) {
	for (int y = 0; y < m_height; y++) {
		for (int x = 0; x < m_width; x++) {
			int index = (y * m_width + x) * 3;
			double value = data[y * m_width + x];
			s_rgb color = getColor(value, scale);
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
s_rgb GenImage::getColor(double value, double supremum) {
	s_rgb color;
	double blue_value =  (value < 0.0) ? ( - value / supremum) : 0.0; //(supremum - value) / (2.0 * supremum);
	double red_value = (value > 0.0) ? (value / supremum) : 0.0; //(supremum + value) / (2.0 * supremum);
	double green_value = 1.0 - (abs(value) / supremum);
	if (blue_value > 1.0 || red_value > 1.0 || green_value > 1.0) {
		blue_value = 0.0;
		red_value = 0.0;
		green_value = 0.0;
	}
	else if (blue_value < 0.0 || red_value < 0.0 || green_value < 0.0) {
		blue_value = 0.0;
		red_value = 0.0;
		green_value = 0.0;
	}
	if (supremum > 0.0) {
		color.c[2] = static_cast<unsigned char>(std::round(255.0 * red_value)); // Red
		color.c[1] = static_cast<unsigned char>(std::round(255.0 * green_value)); // Green
		color.c[0] = static_cast<unsigned char>(std::round(255.0 * blue_value)); // Blue
	}
	else {
		color.c[0] = 0; // Blue
		color.c[1] = 0; // Green
		color.c[2] = 0; // Red
	}
	return color;
}
s_rgb GenImage::getColorRB(double value, double supremum) {
	s_rgb color;
	if (supremum > 0.0) {
		double normValue = value / supremum;
		color.c[2] = static_cast<unsigned char>(std::round(255.0 * normValue)); // Red
		color.c[1] = 0; // Green
		color.c[0] = static_cast<unsigned char>(std::round(255.0 * (1.0 - normValue))); // Blue
	}
	else {
		color.c[0] = 0; // Blue
		color.c[1] = 0; // Green
		color.c[2] = 0; // Red
	}
	return color;
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
bool GenImage::isInBlowGridLine(int x, int y) {
	int half_grid_thickness = m_blow_grid_thickness / 2;
	int x_start = x + half_grid_thickness;
	int y_start = y + half_grid_thickness;
	return (x_start % m_blow_factor < m_blow_grid_thickness) || (y_start % m_blow_factor < m_blow_grid_thickness);
}
void GenImage::drawBlowSingleColorGrid(s_rgb color_foreground, s_rgb color_background) {
	if(m_blow_grid_thickness <= 0 || m_blow_factor <= 1)
		return;
	for (int y = 0; y < m_height; y++) {
		for (int x = 0; x < m_width; x++) {
			int index = (y * m_width + x) * 3;
			bool isGridLine = isInBlowGridLine(x,y);
			s_rgb color = isGridLine ? color_foreground : color_background;
			m_pImageData[index] = color.c[0]; // Red
			m_pImageData[index + 1] = color.c[1]; // Green
			m_pImageData[index + 2] = color.c[2]; // Blue
		}
	}
}
bool GenImage::isInImage(int x, int y) {
	if (x < 0 || y < 0)
		return false;
	if (x >= m_width || y >= m_height)
		return false;
	return true;
}