#include "drawTest.h"

drawTest::drawTest(int grid_width, int grid_height): m_mouse_clicks(0){
	m_pGenImage = new GenImage(grid_width, grid_height);
}
drawTest::~drawTest() {
	if (m_pGenImage != nullptr)
		delete m_pGenImage;
}

void drawTest::drawData(double* data, double dat_sup) {
	m_pGenImage->genScaledImage(data, dat_sup);
}
