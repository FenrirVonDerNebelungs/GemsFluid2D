#include "drawVelocity.h"
drawVelocity::drawVelocity(
	int width_unblown,
	int height_unblown,
	double delta_x,
	double delta_t,
	unsigned char* blow_background_image,
	int blow_factor,
	int blow_grid_thickness,
	int marker_size,
	double display_scale_factor
):m_width_unblown(width_unblown), 
m_height_unblown(height_unblown),  
m_blow_factor(blow_factor),
m_pImage_velocity(nullptr), 
m_pImage_backVelocity(nullptr),
m_pImage_displacement(nullptr),
m_pImage_newVelocity(nullptr),
m_delta_x(delta_x),
m_delta_t(delta_t),
m_marker_size(marker_size),
m_display_scale_factor(display_scale_factor)
{
	if (blow_factor <= 1)
		blow_factor = 1;
	setDefaultColors();/*sets the s_rgb colors*/
	int width = width_unblown * blow_factor;
	int height = height_unblown * blow_factor;
	m_pGenImage = new GenImage(width, height, blow_factor, blow_grid_thickness);
	for (int i = 0; i < 3; i++) {
		m_pImage_3bilinearUx[i] = nullptr;
		m_pImage_3bilinearUy[i] = nullptr;
	}
}
drawVelocity::~drawVelocity() {
	release();
	if (m_pGenImage != nullptr)
		delete m_pGenImage;
}
bool drawVelocity::init(
	const double* Ux_start,
	const double* Uy_start,
	const double* relPos_i,
	const double* relPos_j,
	const double* Ux_new,
	const double* Uy_new,
	const double* bilinear_Ux_approx,
	const double* bilinear_Uy_approx
)
{
	size_t size_of_blown_image_in_data = m_pGenImage->getWidthHeight().width * m_pGenImage->getWidthHeight().height * 3;
	int size_in_pix_of_blown_image = m_pGenImage->getWidthHeight().width * m_pGenImage->getWidthHeight().height;
	m_pImage_velocity = new unsigned char[size_of_blown_image_in_data];
	m_pImage_backVelocity = new unsigned char[size_of_blown_image_in_data];
	m_pImage_displacement = new unsigned char[size_of_blown_image_in_data];
	m_pImage_newVelocity = new unsigned char[size_of_blown_image_in_data];

	draw(Ux_start, Uy_start);
	m_pGenImage->makeCopyOfImageData(m_pImage_velocity);
	draw(Ux_start, Uy_start, false);
	m_pGenImage->makeCopyOfImageData(m_pImage_backVelocity);
	drawDisplacementMarkers(relPos_i, relPos_j);
	m_pGenImage->makeCopyOfImageData(m_pImage_displacement);
	m_end_color.c[2] = 0xff;
	m_end_color.c[1] = 0x0a;
	m_end_color.c[0] = 0xec;
	draw(Ux_new, Uy_new, true, relPos_i, relPos_j, false);
	m_pGenImage->makeCopyOfImageData(m_pImage_newVelocity);

	if (bilinear_Ux_approx != nullptr && bilinear_Uy_approx!=nullptr) {
		for (int i = 0; i < 3; i++) {
			m_pImage_3bilinearUx[i] = new unsigned char[size_of_blown_image_in_data];
			m_pImage_3bilinearUy[i] = new unsigned char[size_of_blown_image_in_data];
		}

		double bilinear_UxUy_max = 0.0;
		for (int i = 0; i < size_in_pix_of_blown_image; i++) {
			if (abs(bilinear_Ux_approx[i]) > bilinear_UxUy_max)
				bilinear_UxUy_max = abs(bilinear_Ux_approx[i]);
			if (abs(bilinear_Uy_approx[i]) > bilinear_UxUy_max)
				bilinear_UxUy_max = abs(bilinear_Uy_approx[i]);
		}
		m_pGenImage->genScaledImage(bilinear_Ux_approx, bilinear_UxUy_max);
		m_pGenImage->makeCopyOfImageData(m_pImage_3bilinearUx[0]);
		drawDisplacementMarkers(relPos_i, relPos_j);
		m_pGenImage->makeCopyOfImageData(m_pImage_3bilinearUx[1]);
		drawDisplacementMarkers(relPos_i, relPos_j, Ux_new, bilinear_UxUy_max);
		m_pGenImage->makeCopyOfImageData(m_pImage_3bilinearUx[2]);

		m_pGenImage->genScaledImage(bilinear_Uy_approx, bilinear_UxUy_max);
		m_pGenImage->makeCopyOfImageData(m_pImage_3bilinearUy[0]);
		drawDisplacementMarkers(relPos_i, relPos_j);
		m_pGenImage->makeCopyOfImageData(m_pImage_3bilinearUy[1]);
		drawDisplacementMarkers(relPos_i, relPos_j, Uy_new, bilinear_UxUy_max);
		m_pGenImage->makeCopyOfImageData(m_pImage_3bilinearUy[2]);
	}
	return true;
}
void drawVelocity::release() {
	for (int i = 0; i < 3; i++) {
		if (m_pImage_3bilinearUx[i] != nullptr) {
			delete[] m_pImage_3bilinearUx[i];
			m_pImage_3bilinearUx[i] = nullptr;
		}
		if (m_pImage_3bilinearUy[i] != nullptr) {
			delete[] m_pImage_3bilinearUy[i];
			m_pImage_3bilinearUy[i] = nullptr;
		}
	}
	if (m_pImage_velocity != nullptr) {
		delete []m_pImage_velocity;
		m_pImage_velocity = nullptr;
	}
	if (m_pImage_backVelocity != nullptr) {
		delete []m_pImage_backVelocity;
		m_pImage_backVelocity = nullptr;
	}
	if (m_pImage_displacement != nullptr) {
		delete[] m_pImage_displacement;
		m_pImage_displacement = nullptr;
	}
	if (m_pImage_newVelocity != nullptr) {
		delete[] m_pImage_newVelocity;
		m_pImage_newVelocity = nullptr;
	}
}
void drawVelocity::draw(const double* Ux, const double* Uy, bool draw_forward, const double* i_offset, const double* j_offset, bool draw_background)
{
	double sign = draw_forward ? 1.0 : -1.0;
	if(draw_background)
		m_pGenImage->drawBlowSingleColorGrid(m_grid_color, m_background_color);
	for (int j = 0; j < m_height_unblown; j++) {
		for (int i = 0; i < m_width_unblown; i++) {
			int index = j * m_width_unblown + i;
			/*get velocity in original grid screen scale*/
			double di = Ux[index]*m_delta_t/m_delta_x;
			double dj = Uy[index]*m_delta_t/m_delta_x;
			/*reverse sign if needed*/
			di *= sign;
			dj *= sign;
			/*get start location in blown screen coord from unblown cell center coord */
			int x_start, y_start;
			if(!getBlowScreenCoord_from_coord(x_start, y_start, i, j)) {
				continue;/*if start location is out of screen, skip*/
			}
			/*get end location after displacement in screen coord of di, dj*/
			int x_end, y_end;
			if (!getBlowScreenCoord_from_coordDisplacement(x_end, y_end, i, j, di, dj)) {
				continue;/*if end location is out of screen, skip*/
			}
			/* if there are offsets relocate start in blown screen coordinates */
			if(i_offset!=nullptr) {
				double di_offset = i_offset[index];
				x_start += static_cast<int>(std::round(di_offset * m_pGenImage->getBlowFactor()));
				x_end += static_cast<int>(std::round(di_offset * m_pGenImage->getBlowFactor()));
			}
			if(j_offset!=nullptr) {
				double dj_offset = j_offset[index];
				y_start += static_cast<int>(std::round(dj_offset * m_pGenImage->getBlowFactor()));
				y_end += static_cast<int>(std::round(dj_offset * m_pGenImage->getBlowFactor()));
			}

			if(draw_forward)
				m_pGenImage->drawGradientLine(x_start, y_start, x_end, y_end, m_start_color, m_end_color);
			else
				m_pGenImage->drawGradientLine(x_start, y_start, x_end, y_end, m_end_color, m_start_color);
		}
	}
}
void drawVelocity::drawBackwards(const double* Ux, const double* Uy) {
	draw(Ux, Uy, false);
}

void drawVelocity::drawDisplacementMarkers(const double* di, const double* dj, const double* marker_col, double supremum) {
	for (int j = 0; j < m_height_unblown; j++) {
		for (int i = 0; i < m_width_unblown; i++) {
			int index = j * m_width_unblown + i;
			s_rgb color = m_marker_color;
			if(marker_col!=nullptr) {
				color = getColorFromMarkerCol(marker_col[index], supremum);
			}
			int x_blow, y_blow;
			if(getBlowScreenCoord_from_coordDisplacement(x_blow, y_blow, i, j, di[index], dj[index])) {
				m_pGenImage->drawSolidBox(x_blow, y_blow, m_marker_size, color);
			}
		}
	}
}
bool drawVelocity::getBlowScreenCoord_from_coord(int& i_blow, int& j_blow, int i, int j) {
	double blow_factor = (double)m_pGenImage->getBlowFactor();
	double half_square = blow_factor / 2.0;
	double x_center = (double)i;
	double y_center = (double)j;
	double x_center_blow = x_center * blow_factor + half_square;
	double y_center_blow = y_center * blow_factor + half_square;
	i_blow = static_cast<int>(std::round(x_center_blow));
	j_blow = static_cast<int>(std::round(y_center_blow));
	return (i_blow >= 0 && i_blow < m_pGenImage->getWidthHeight().width && j_blow >= 0 && j_blow < m_pGenImage->getWidthHeight().height);
}
bool drawVelocity::getBlowScreenCoord_from_coordDisplacement(int& i_blow, int& j_blow, int i, int j, double di, double dj) {
	double blow_factor = (double)m_pGenImage->getBlowFactor();
	double half_square = blow_factor / 2.0;
	double x_center = (double)i;
	double y_center = (double)j;
	double x_center_blow = x_center * blow_factor + half_square;
	double y_center_blow = y_center * blow_factor + half_square;
	/*scale the dist vector to blowup scale*/
	double di_blow = di * m_pGenImage->getBlowFactor()*m_display_scale_factor;
	double dj_blow = dj * m_pGenImage->getBlowFactor()*m_display_scale_factor;
	/*debug*/
	if (di > 0.0)
		int testttt = 0;
	/*     */
	double x_blow = x_center_blow + di_blow;
	double y_blow = y_center_blow + dj_blow;
	i_blow = static_cast<int>(std::round(x_blow));
	j_blow = static_cast<int>(std::round(y_blow));
	return (i_blow >= 0 && i_blow < m_pGenImage->getWidthHeight().width && j_blow >= 0 && j_blow < m_pGenImage->getWidthHeight().height);
}
s_rgb drawVelocity::getColorFromMarkerCol(double marker_col, double supremum) {
	return m_pGenImage->getColor(marker_col, supremum);
}

void drawVelocity::setDefaultColors() {
	for (int i = 0; i < 3; i++)
		m_grid_color.c[i] = 0x05;
	m_background_color.c[2] = 0xcc;
	m_background_color.c[1] = 0xe6;
	m_background_color.c[0] = 0xFF;
	m_start_color.c[2] = 0xcc;
	m_start_color.c[1] = 0xe6;
	m_start_color.c[0] = 0xFF;
	m_end_color.c[2] = 0xff;
	m_end_color.c[1] = 0x1a;
	m_end_color.c[0] = 0x1a;
	m_marker_color.c[2] = 0xff;
	m_marker_color.c[1] = 0x1a;
	m_marker_color.c[0] = 0xff;

}