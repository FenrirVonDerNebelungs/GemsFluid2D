#pragma once
#ifndef DRAWTEST_H
#define DRAWTEST_H
#ifndef GENIMAGE_H
#include "GenImage.h"
#endif
#include <string>

class Test;

class drawTest {
public:
	drawTest(int grid_width, int grid_height);
	~drawTest();
	const WCHAR* GetMessage() const { return (message_str.length()>1) ? message_str.c_str() : nullptr; }
	void drawData(double* data, double dat_sup);

	GenImage* getGenImage() { return m_pGenImage; }
	friend class Test;
protected:
	GenImage* m_pGenImage;
	int m_mouse_clicks;
	std::wstring message_str;

	void setMessage(std::wstring mes) { message_str=mes; }
};
#endif
