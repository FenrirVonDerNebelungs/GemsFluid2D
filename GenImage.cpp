#include "GenImage.h"

GenImage::GenImage(int width, int height)
{
	m_width = width;
	m_height = height;
	int nBitsPerPixel = 24;
	//int nBytesPerPixel = nBitsPerPixel / 8;
	//int nBytesPerLine = ((m_width * nBytesPerPixel + 3) / 4) * 4; 
	//int nImageSize = nBytesPerLine * m_height;
	m_pbmi = reinterpret_cast<BITMAPINFO*>( GlobalAlloc(GPTR, sizeof(BITMAPINFOHEADER)) ); // +sizeof(RGBQUAD) * colorsUsed);  GPTR is GMEM_FIXED | GMEM_ZEROINIT which means bytes are in fixed memory referenced by a pointer and mem is zeroed
	m_pbmi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	m_pbmi->bmiHeader.biWidth = m_width;
	m_pbmi->bmiHeader.biHeight = m_height;//-m_height;  autcomplete made this negative and added the comment "top-down image" see what happens with this positive
	m_pbmi->bmiHeader.biPlanes = 1;
	m_pbmi->bmiHeader.biBitCount = nBitsPerPixel;
	m_pbmi->bmiHeader.biCompression = BI_RGB;/* no compression*/
	m_pbmi->bmiHeader.biSizeImage = 0;/* can be 0 for BI_RGB (not compressed) bitmaps*/
	m_pbmi->bmiHeader.biXPelsPerMeter = 0;
	m_pbmi->bmiHeader.biYPelsPerMeter = 0;
	m_pbmi->bmiHeader.biClrUsed = 0;/* we are using all 16 million colors; use the maximum number of colors for the bitdepth 2^32 */
	m_pbmi->bmiHeader.biClrImportant = 0;

	m_pImageData = new unsigned char[m_width * m_height * nBitsPerPixel / 8];
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
	if(data==nullptr)
		return;
	double supremum = getDataSupremum(data);
	for (int y = 0; y < m_height; y++) {
		for (int x = 0; x < m_width; x++) {
			int index = (y * m_width + x) * 3;
			double value = data[y * m_width + x];
			s_rgb color = getColor(value, supremum);
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
		if (data[i] > supremum) {
			supremum = data[i];
		}
	}
	return supremum;
}
s_rgb GenImage::getColor(double value, double supremum) {
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