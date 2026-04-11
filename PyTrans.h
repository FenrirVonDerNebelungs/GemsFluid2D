#pragma once
#ifndef PYTRANS_H
#define PYTRANS_H

#ifndef BASE_H
#include "Base.h"
#endif
#include<stdint.h>

#define PYTRANS_NUMHEADER_FIELDS 11
#define PYTRANS_FILE_HEADER_LEN 7
#define PYTRANS_CACHE_HEADER_LEN 7

namespace n_PyTrans {
	/*These are the IDs for each object type used to identify the object transmitted, the int32_t is sent at the start of the char string
	to read this these codes will have to be replicated in the python code*/
	/*1st 32 bits*/
	const int32_t grid_code = 0x01;
	const int32_t img_code = 0x02;
	/*second 32 bits*/
	const int32_t U_code = 0x01;
	const int32_t Force_code = 0x02;
	const int32_t relPos_code = 0x03;
	const int32_t W_code = 0x04;
	const int32_t DivW_code = 0x05;
	const int32_t P_code = 0x06;/*pressure*/
	const int32_t gradP_code = 0x07;
	const int32_t jacobi_frame_code = 0x08;
	const int32_t DivU_code = 0x09;
	/* 3rd 32 bits */
	const int32_t X_code = 0x01;
	const int32_t Y_code = 0x02;
	const int32_t Scalar_code = 0x00;/*scalar*/
	const int32_t X_exp_code = 0x03;/*enlarged */
	const int32_t Y_exp_code = 0x04;
	/*4th 32 bits step code */
	const int32_t mid_frame_code = 0x00;
	const int32_t start_frame_code = 0x01;
	const int32_t after_advection_code = 0x02;
	const int32_t after_force_code = 0x03;
	const int32_t end_frame_code = 0x04;
	/*5th 32 bit code 
		number of loops processed */
	/*6th 32 bits cnt code
		jacobi frame count */
	/* 7th 32 bits cnt code
		grid width */
	/* 8th 32 bits cnt code
		grid height */
	/* 9th 32 bits cnt code*
		grid_expansion
	*/
	/* 10th 32 bits code
		max value in grid or image, used for normalization in python 
		*/
	/* 11th 32 bits code
		min value in grid or image, used for normalization in python
		*/ 
	const int obj_header_len = PYTRANS_NUMHEADER_FIELDS;

	/*file header*/
	/*1st 32 bits
	* len of file header in bytes (chars)
	*/
	/* 2nd 32 bits
	* len of each cache header in bytes (chars)
	*/
	/* 3rd 32 bits 
	* len of each stream header
	*/
	/*4th 32 bits
	* number of caches written
	*/
	/* 5th 32 bits
	* max number of headers per cache
	*/
	/* 6th 32 bits
	* max number of stream headers per cache
	*/
	/* 7th 32 bits
	* max number of image headers per cache
	*/
	const int file_header_len = PYTRANS_FILE_HEADER_LEN;

	/*1st len of cache header in bytes*/
	/*2nd len of each stream header in bytes */
	/*3rd number of headers per cache*/
	/*4th number of stream header per cache*/
	/*5th number of image header per cache*/ 
	/*6th data type for stream
	* 0x01 is float 32
	* 0x02 is double 64
	*/
	/*7th cache size in bytes, after header*/
	const int32_t data_type_float_code = 0x01;
	const int32_t data_type_double_code = 0x02;

	const int cache_header_len = PYTRANS_CACHE_HEADER_LEN;

	/* conversion functions used so that int and float data can be packed into a char array which is then either sent or
	   recieved by python, it is assumed that 'float' is a 32 bit float that follows the IEEE 754 standard, this is supposed to be the default standard for windows*/
	void charTochar4(int char_loc, const char chs[], char ch4[]);/*starts at char_loc and takes 4 values from chs and puts them in ch4 which should have a length of 4*/
	void char4Tochar(const char ch4[], int char_loc, char chs[]);/*takes the 4 values from ch4 and inserts them at char_loc into chs*/
	void charTochar8(int char_loc, const char chs[], char ch8[]);
	void char8Tochar(const char ch8[], int char_loc, char chs[]);
	void char3Tochar4(const char ch3[], char ch4[]);
	void char4Tochar3(const char ch4[], char ch3[]);
	int32_t Fti(float F);/*direct bit cast to an int32_t, used for transmission, the code at the other end should directly reverse the cast to put it back into a 32bit IEEE 754 float*/
	int64_t Dti(double D);
	float Itf(int32_t I);/*direct bit cast from int32_t to a float*/
	double Itd(int64_t I);
	bool Itc(const int32_t I, char ch[]);/*puts the int32_t into the char array chars by default should be 8 bits, or 1 byte apiece, array must be 4 chars long */
	bool Ltc(const int64_t I, char ch[]);
	bool Cti(const char ch[], int32_t& I);/*puts the char array into the int32, array must be 4 chars long*/
	bool Ctl(const char ch[], int64_t& I);
	int sizeIFromChar(int len_char); 
	int sizeLFromChar(int len_char);
	int sizeCharFromI(int len_I);
	int sizeCharFromL(int len_I);
	bool streamItoC(int len, const int32_t Is[], char chs[], int char_offset = 0);/*the ch array must have 4 times the length of the int array*/
	bool streamLtoC(int len, const int64_t Is[], char chs[], int char_offset = 0);
	bool streamCtoI(int len, const char chs[], int32_t Is[], int char_offset = 0);/* the len is the len of the ch array INCLUDING the offset */
	bool streamCtoL(int len, const char chs[], int64_t Is[], int char_offset = 0);
	bool streamDtoC(int len, const double Fs[], char chs[], int char_offset = 0);
	bool streamCtoD(int len, const char chs[], double Fs[], int char_offset = 0);
}

class PyTrans {
public:
	PyTrans();
	~PyTrans();

	bool init(const char* filename=nullptr, int total_stream_len_in_doubles = 0, int num_headers=1, int max_img_size_in_pix=0, int num_images=1);
	bool releaseAndWrite();
	void release();

	bool cacheGrid(
		const double grid_vals[],
		int grid_width,
		int grid_height,
		int number_of_loops=0,
		int32_t label_code = 0,
		int32_t axis_code=0,
		int32_t start_end_code = 0,
		int jacobi_frame=0
	);
	bool cacheImage(
		const unsigned char img_vals[], /*assumes 3 char per pix input stream size is width*height*3 */
		int img_width,
		int img_height,
		int number_of_loops = 0,
		int32_t label_code = 0,
		int32_t axis_code = 0,
		int32_t start_end_code = 0,
		float max = 0.f,
		float min = 0.f,
		int expansion = 1
	);
	bool cacheDStream(const double d_vals[], int len);
	void resetCache();

	bool resetAndWrite();
private:
	int m_file_header_len;/* len file header in chars (bytes)*/
	int m_cache_header_len;/*len of each cache header in chars (bytes)*/
	int m_header_len;/*stream header length in chars, not int32_t's */
	int m_num_caches_written;
	int m_num_grid_headers_cached;/*per cache */
	int m_num_image_headers_cached;/*per cache */
	int m_max_grid_headers_cached;
	int m_max_image_headers_cached;
	char* m_cache_header;
	char* m_culmative_stream;
	int   m_culmative_stream_len;
	int   m_culmative_stream_len_max;
	char* m_image_stream;
	int   m_image_stream_len;
	int   m_image_stream_len_max;
	char  m_filename[30];


	bool writeFileHeader();
	bool writeBinary();

	bool cacheHeader();

	bool sendGrid(
		char ch_stream[],
		const double grid_vals[],
		int grid_width,
		int grid_height,
		int32_t label_code = 0,
		int32_t axis_code = 0,
		int32_t start_end_code = 0,
		int number_of_loops = 0,
		int jacobi_frame = 0,
		int ch_stream_offset=0,
		int ch_stream_len = -1);
	int getGridStreamLen(int grid_width, int grid_height);/*assumes double for grid includes header len*/

	bool sendImage(
		char ch_stream[],
		const unsigned char img_vals[],
		int img_width,
		int img_height,
		int32_t label_code = 0,
		int32_t axis_code = 0,
		int32_t start_end_code = 0,
		int number_of_loops = 0,
		float max = 0.f,
		float min = 0.f,
		int expansion = 1,
		int ch_stream_offset = 0,
		int ch_stream_len = -1
	);

	bool sendCacheHeader(
		char ch_stream[],
		int number_of_stream_headers,
		int number_of_image_headers,
		int data_type_code,
		int cache_size_in_bytes,
		int ch_stream_offset=0,
		int ch_stream_len=-1
	);

	bool sendHeader(
		char ch_stream[],
		int32_t type_code,
		int32_t label_code,
		int32_t axis_code,
		int32_t start_end_code,
		int number_of_loops,
		int jacobi_frame,
		int width,
		int height,
		int expansion,
		float max,
		float min,
		int ch_stream_offset=0,
		int ch_stream_len=-1
	);
	
};
#endif