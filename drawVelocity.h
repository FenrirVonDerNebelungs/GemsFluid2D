#pragma once
#ifndef DRAWVELOCITY_H
#define DRAWVELOCITY_H

#ifndef GENIMAGE_H
#include "GenImage.h"
#endif

class drawVelocity
{
public:
	drawVelocity(
		int width_unblown, 
		int height_unblown, 
		double delta_x, 
		double delta_t, 
		unsigned char* blow_background_image = nullptr, 
		int blow_factor = 6, 
		int blow_grid_thickness = 1,
		int marker_size=4,
		double display_scale_factor=1
	);
	~drawVelocity();
	bool init(
		const double* Ux_start,
		const double* Uy_start,
		const double* relPos_i,
		const double* relPos_j,
		const double* Ux_new,
		const double* Uy_new,
		const double* bilinear_Ux_aprox = nullptr,
		const double* bilinear_Uy_aprox = nullptr
	);
	void release();

	void draw(
		const double* Ux /*velocity in real space*/,
		const double* Uy /*velocity in real space*/,
		bool draw_forward = true /*foward starts at offset position in center of sell and draws vector from there, backwards draws to */,
		const double* i_offset = nullptr /* offset distance for x in unblown screen coordinates */,
		const double* j_offset = nullptr /* offset distance for y in unblown screen coordinates */,
		bool draw_background = true);/* if true draws grid square background, if false just draws velocity*/
	void drawBackwards(const double* Ux, const double* Uy);/*draws velocity vectors so that they end in the cell that has
																					   their velocity */
	void drawDisplacementMarkers(const double* di, const double* dj, const double* marker_col=nullptr, double supremum=1.0);

	bool setImage_to_velocity() { return m_pGenImage->copyImageDataFromBuffer(m_pImage_velocity); }
	bool setImage_to_backVelocity() { return m_pGenImage->copyImageDataFromBuffer(m_pImage_backVelocity); }
	bool setImage_to_displacement() { return m_pGenImage->copyImageDataFromBuffer(m_pImage_displacement); }
	bool setImage_to_newVelocity() { return m_pGenImage->copyImageDataFromBuffer(m_pImage_newVelocity); }
	bool setImage_to_bilinear_Ux(int index) { return m_pGenImage->copyImageDataFromBuffer(m_pImage_3bilinearUx[index]); }
	bool setImage_to_bilinear_Uy(int index) { return m_pGenImage->copyImageDataFromBuffer(m_pImage_3bilinearUy[index]); }

	GenImage* getGenImage() { return m_pGenImage; }
	int getBlowFactor() { return m_blow_factor; }
private:
	int m_width_unblown;
	int m_height_unblown;
	int m_blow_factor;
	GenImage* m_pGenImage;
	unsigned char* m_pImage_velocity;
	unsigned char* m_pImage_backVelocity;
	unsigned char* m_pImage_displacement;/*displacement estimate set on top of back velocity*/
	unsigned char* m_pImage_newVelocity;/*new velocity estimate set on top of displacement marker set on top of prev foward velocity*/
	unsigned char* m_pImage_3bilinearUx[3];/* the first image is the bilinear interpolation of Ux in the blown screen coordinates,
	                                          the second is the image of the bilinear interpolation with the markers at the backpropogated velocity tail starts
											  the third is the same image but with the markers colored by the Ux computed by Advection
	                                          if Advection works properly the colors of the final markers should match the colors computed by the bilinear inerpolation*/
	unsigned char* m_pImage_3bilinearUy[3];/* same as above but with Uy*/
	double m_delta_x;
	double m_delta_t;
	s_rgb m_grid_color;
	s_rgb m_background_color;
	s_rgb m_start_color;
	s_rgb m_end_color;
	s_rgb m_marker_color;
	int   m_marker_size;
	double m_display_scale_factor;

	bool getBlowScreenCoord_from_coord(int& i_blow, int& j_blow, int i, int j);
	bool getBlowScreenCoord_from_coordDisplacement(int& i_blow, int& j_blow, int i, int j, double di, double dj);

	s_rgb getColorFromMarkerCol(double marker_col, double supremum=1.0);

	void setDefaultColors();
};

#endif