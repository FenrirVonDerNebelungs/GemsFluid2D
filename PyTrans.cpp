#include "PyTrans.h"

#include <fstream>
//#include <bit>

void n_PyTrans::charTochar4(int char_loc, const char chs[], char ch4[]) {
	int chs_i = char_loc;
	for (int i = 0; i < 4; i++) {
		ch4[i] = chs[chs_i];
		chs_i++;
	}
}
void n_PyTrans::char4Tochar(const char ch4[], int char_loc, char chs[]) {
	int chs_i = char_loc;
	for (int i = 0; i < 4; i++) {
		chs[chs_i] = ch4[i];
		chs_i++;
	}
}
void n_PyTrans::charTochar8(int char_loc, const char chs[], char ch8[]) {
	int chs_i = char_loc;
	for (int i = 0; i < 8; i++) {
		ch8[i] = chs[chs_i];
		chs_i++;
	}
}
void n_PyTrans::char8Tochar(const char ch8[], int char_loc, char chs[]) {
	int chs_i = char_loc;
	for (int i = 0; i < 8; i++) {
		chs[chs_i] = ch8[i];
		chs_i++;
	}
}
void n_PyTrans::char3Tochar4(const char ch3[], char ch4[]) {
	for (int i = 0; i < 3; i++)
		ch4[i] = ch3[i];
	ch4[3] = 0x00;
}
void n_PyTrans::char4Tochar3(const char ch4[], char ch3[]) {
	for (int i = 0; i < 3; i++)
		ch3[i] = ch4[i];
}
int32_t n_PyTrans::Fti(float F) {
	if (sizeof(F) != 4)
		return std::numeric_limits<int32_t>::min();
	int32_t convval = std::_Bit_cast<int32_t, float>(F);//bit_cast<int32_t>(F);
	return convval;
}
int64_t n_PyTrans::Dti(double D) {
	if (sizeof(D) != 8)
		return std::numeric_limits<int64_t>::min();
	int64_t convval = std::_Bit_cast<int64_t, double>(D);
	return convval;
}
float n_PyTrans::Itf(int32_t I) {
	float convval = std::_Bit_cast<float, int32_t>(I);//bit_cast<float>(I);
	if (sizeof(convval) != 4)
		return std::numeric_limits<float>::quiet_NaN();
	return convval;
}
double n_PyTrans::Itd(int64_t I) {
	double convval = std::_Bit_cast<double, int64_t>(I);
	if (sizeof(convval) != 8)
		return std::numeric_limits<double>::quiet_NaN();
	return convval;
}
bool n_PyTrans::Itc(const int32_t I, char ch[]) {
	/* 00 00 00 00 */
	int32_t iext = I >> 24;
	ch[3] = static_cast<char>(iext);
	iext = I >> 16;
	int32_t mI = iext & 0xFF;
	ch[2] = static_cast<char>(mI);
	iext = I >> 8;
	mI = iext & 0xFF;
	ch[1] = static_cast<char>(mI);
	mI = I & 0xFF;
	ch[0] = static_cast<char>(mI);
	return true;
}
bool n_PyTrans::Ltc(const int64_t I, char ch[]) {
	int64_t hi = I >> 32;
	int32_t hi_I = static_cast<int32_t>(hi);
	int64_t lo = I & 0xFFFFFFFF;
	int32_t lo_I = static_cast<int32_t>(lo);
	char hi_c[4];
	char lo_c[4];
	if (!Itc(hi_I, hi_c))
		return false;
	if (!Itc(lo_I, lo_c))
		return false;
	for(int i=0; i<4; i++){
		ch[i] = lo_c[i];
		ch[i + 4] = hi_c[i];
	}
	return true;
}

bool n_PyTrans::Cti(const char ch[], int32_t& I) {
	int32_t iext = static_cast<int32_t>(ch[3]);
	I = iext << 24;
	iext = static_cast<int32_t>(ch[2]);
	int32_t ish = iext << 16;
	I = I | ish;
	iext = static_cast<int32_t>(ch[1]);
	ish = iext << 8;
	I = I | ish;
	iext = static_cast<int32_t>(ch[0]);
	I = I | iext;
	return true;
}
bool n_PyTrans::Ctl(const char ch[], int64_t& I) {
	char hi_c[4];
	char lo_c[4];
	for (int i = 0; i < 4; i++) {
		lo_c[i] = ch[i];
		hi_c[i] = ch[i + 4];
	}
	int32_t lo_I, hi_I;
	if (!Cti(lo_c, lo_I))
		return false;
	if (!Cti(hi_c, hi_I))
		return false;
	int64_t lo = static_cast<int64_t>(lo_I);
	int64_t hi = static_cast<int64_t>(hi_I);
	I = hi << 32;
	I = I | lo;
	return true;
}
int n_PyTrans::sizeIFromChar(int len_char) { return len_char / 4; }
int n_PyTrans::sizeLFromChar(int len_char) { return len_char / 8; }
int n_PyTrans::sizeCharFromI(int len_I) { return len_I * 4; }
int n_PyTrans::sizeCharFromL(int len_I) { return len_I * 8; }
bool n_PyTrans::streamItoC(int len, const int32_t Is[], char chs[], int char_offset) {
	char ch4[4];
	int chs_i = char_offset;
	for (int i = 0; i < len; i++) {
		Itc(Is[i], ch4);
		char4Tochar(ch4, chs_i, chs);
		chs_i += 4;
	}
	return true;
}
bool n_PyTrans::streamLtoC(int len, const int64_t Is[], char chs[], int char_offset) {
	char ch8[8];
	int chs_i = char_offset;
	for (int i = 0; i < len; i++) {
		Ltc(Is[i], ch8);
		char8Tochar(ch8, chs_i, chs);
		chs_i += 8;
	}
	return true;
}
bool n_PyTrans::streamDtoC(int len, const double Fs[], char chs[], int char_offset) {
	char ch8[8];
	int chs_i = char_offset;
	for (int i = 0; i < len; i++) {
		int64_t d_I = Dti(Fs[i]);
		Ltc(d_I, ch8);
		char8Tochar(ch8, chs_i, chs);
		chs_i += 8;
	}
	return true;
}
bool n_PyTrans::streamCtoI(int len, const char chs[], int32_t Is[], int char_offset) {
	char ch4[4];
	int Is_i = 0;
	int32_t Int_to_unpack = 0;
	for (int chs_i = 0; chs_i < len; chs_i += 4) {
		charTochar4(chs_i+char_offset, chs, ch4);
		Cti(ch4, Int_to_unpack);
		Is[Is_i] = Int_to_unpack;
		Is_i++;
	}
	return true;
}
bool n_PyTrans::streamCtoL(int len, const char chs[], int64_t Is[], int char_offset) {
	char ch8[8];
	int Is_i = 0;
	int64_t Int_to_unpack = 0;
	for (int chs_i = 0; chs_i < len; chs_i += 8) {
		charTochar8(chs_i+char_offset, chs, ch8);
		Ctl(ch8, Int_to_unpack);
		Is[Is_i] = Int_to_unpack;
		Is_i++;
	}
	return true;
}
bool n_PyTrans::streamCtoD(int len, const char chs[], double Fs[], int char_offset) {
	char ch8[8];
	int Fs_i = 0;
	int64_t I_to_unpack = 0;
	for (int chs_i = 0; chs_i < len; chs_i += 8) {
		charTochar8(chs_i+char_offset, chs, ch8);
		Ctl(ch8, I_to_unpack);
		double d_to_unpack = Itd(I_to_unpack);
		Fs[Fs_i] = d_to_unpack;
		Fs_i++;
	}
	return true;
}
PyTrans::PyTrans() :
	m_num_caches_written(0),
	m_num_grid_headers_cached(0),
	m_num_image_headers_cached(0),
	m_max_grid_headers_cached(0),
	m_max_image_headers_cached(0),
	m_cache_header(nullptr),
	m_culmative_stream(nullptr), 
	m_culmative_stream_len(0), 
	m_culmative_stream_len_max(0),
	m_image_stream(nullptr), 
	m_image_stream_len(0), 
	m_image_stream_len_max(0) 
{
	m_file_header_len = sizeof(int32_t) * n_PyTrans::file_header_len;	/*sizeof returns size in bytes 1 char = 1 byte*/
	m_header_len = n_PyTrans::sizeCharFromI(n_PyTrans::obj_header_len);
	m_cache_header_len = n_PyTrans::sizeCharFromI(n_PyTrans::cache_header_len);
	for (int i = 0; i < 30; i++)
		m_filename[i] = '\0';
}
PyTrans::~PyTrans() {
	;
}
bool PyTrans::init(const char* filename, int total_stream_len_in_doubles, int num_headers, int max_img_size_in_pix, int num_images) {
	bool retval = true;
	if (filename != nullptr)
		strcpy_s(m_filename, filename);
	else
		retval = false;
	int term_loc = 0;
	if(filename!=nullptr)
		term_loc = (int)std::strlen(filename);
	m_filename[term_loc] = '\0';
	m_num_caches_written = 0;
	m_num_grid_headers_cached = 0;
	m_num_image_headers_cached = 0;
	if(m_cache_header_len>0){
		m_cache_header = new char[m_cache_header_len];
		if (m_cache_header == nullptr)
			return false;
		for(int ii=0; ii<m_cache_header_len; ii++)
			m_cache_header[ii] = '\0';
	}
	if (total_stream_len_in_doubles > 0) {
		m_culmative_stream_len_max = n_PyTrans::sizeCharFromL(total_stream_len_in_doubles);
		m_culmative_stream_len_max += (num_headers * m_header_len);
		m_culmative_stream = new char[m_culmative_stream_len_max];
		if (m_culmative_stream == nullptr) {
			release();
			return false;
		}
	}
	if (max_img_size_in_pix > 0) {
		m_image_stream_len_max = 3*sizeof(char) * max_img_size_in_pix;
		m_image_stream_len_max += (num_images * m_header_len);
		m_image_stream = new char[m_image_stream_len_max];
		if (m_image_stream == nullptr) {
			release();
			return false;
		}
	}
	if (total_stream_len_in_doubles <= 0 && max_img_size_in_pix <= 0)
		retval = false;
	if(retval)
		retval = writeFileHeader();
	return retval;
}

void PyTrans::release() {
	m_num_caches_written = 0;
	m_num_grid_headers_cached = 0;
	m_num_image_headers_cached = 0;
	m_max_grid_headers_cached = 0;
	m_max_image_headers_cached = 0;
	if(m_cache_header!=nullptr)
		delete[] m_cache_header;
	m_cache_header = nullptr;
	if(m_image_stream!=nullptr)
		delete[] m_image_stream;
	m_image_stream = nullptr;
	m_image_stream_len_max = 0;
	m_image_stream_len = 0;
	if (m_culmative_stream != nullptr)
		delete[] m_culmative_stream;
	m_culmative_stream = nullptr;
	m_culmative_stream_len_max = 0;
	m_culmative_stream_len = 0;
	m_filename[0] = '\0';
}
bool PyTrans::cacheGrid(
	const double grid_vals[],
	int grid_width,
	int grid_height,
	int number_of_loops,
	int32_t label_code,
	int32_t axis_code,
	int32_t start_end_code,
	int jacobi_frame
) {
	int new_cache_len = getGridStreamLen(grid_width, grid_height);
	if ((new_cache_len + m_culmative_stream_len) >= m_culmative_stream_len_max)
		return false;
	bool retVal = sendGrid(
		m_culmative_stream, 
		grid_vals, 
		grid_width, 
		grid_height, 
		label_code, 
		axis_code, 
		start_end_code, 
		number_of_loops, 
		jacobi_frame, 
		m_culmative_stream_len);
	m_culmative_stream_len += new_cache_len;
	m_num_grid_headers_cached++;
	return retVal;
}
bool PyTrans::cacheImage(
	const unsigned char img_vals[],
	int img_width,
	int img_height,
	int number_of_loops,
	int32_t label_code,
	int32_t axis_code,
	int32_t start_end_code,
	float max,
	float min,
	int expansion
) {
	int image_size = img_width * img_height * 3;
	int new_cache_len = image_size + m_header_len;
	if ((new_cache_len + m_image_stream_len) >= m_image_stream_len_max)
		return false;
	bool retVal = sendImage(
		m_image_stream, 
		img_vals, 
		img_width, 
		img_height, 
		label_code, 
		axis_code, 
		start_end_code, 
		number_of_loops, 
		max, 
		min, 
		expansion, 
		m_image_stream_len
	);
	m_image_stream_len += new_cache_len;
	m_num_image_headers_cached++;
	return retVal;
}
bool PyTrans::cacheDStream(const double d_vals[], int len) {
	if ((len + m_culmative_stream_len) >= m_culmative_stream_len_max)
		return false;
	bool ret_val = n_PyTrans::streamDtoC(len, d_vals, m_culmative_stream, m_culmative_stream_len);
	if (ret_val)
		m_culmative_stream_len += n_PyTrans::sizeCharFromL(len);
	return ret_val;
}
void PyTrans::resetCache() {
	m_culmative_stream_len = 0;
	if(m_num_grid_headers_cached>m_max_grid_headers_cached)
		m_max_grid_headers_cached = m_num_grid_headers_cached;
	if (m_num_image_headers_cached > m_max_image_headers_cached)
		m_max_image_headers_cached = m_num_image_headers_cached;
	m_num_grid_headers_cached = 0;
	m_num_image_headers_cached = 0;
}
bool PyTrans::resetAndWrite() {
	cacheHeader();
	bool good_write = writeBinary();
	resetCache();
	return good_write;
}
bool PyTrans::writeFileHeader() {
	if (m_filename[0] == '\0' )
		return false;
	std::ofstream outF;
	outF.open(m_filename, std::ios::out | std::ios::binary);
	if (!outF.is_open())
		return false;
	char* file_header_ch = new char[m_file_header_len];
	int32_t file_header_I[n_PyTrans::file_header_len];
	file_header_I[0] = m_file_header_len;
	file_header_I[1] = m_cache_header_len;
	file_header_I[2] = m_header_len;
	file_header_I[3] = 0;// m_num_caches_written;
	file_header_I[4] = 0;// (m_max_grid_headers_cached + m_max_image_headers_cached);
	file_header_I[5] = 0; //m_max_grid_headers_cached;
	file_header_I[6] = 0; //m_max_image_headers_cached;
	bool retval = n_PyTrans::streamItoC(n_PyTrans::file_header_len, file_header_I, file_header_ch);
	outF.seekp(0, std::ios::beg);
	outF.write(file_header_ch, m_file_header_len * sizeof(char));
	outF.close();
	delete[] file_header_ch;
	return true;
}
bool PyTrans::writeBinary() {
	if (m_filename[0] == '\0' || (m_culmative_stream_len<=0 && m_image_stream_len<=0) )
		return false;
	std::ofstream outF;
	outF.open(m_filename, std::ios::app | std::ios::binary);
	if (!outF.is_open()) 
		outF.open(m_filename, std::ios::out | std::ios::binary);
	if (!outF.is_open()) {
		return false;
	}
	if(m_cache_header_len>0)
		outF.write(m_cache_header, m_cache_header_len * sizeof(char));
	if (m_culmative_stream_len > 0)
		outF.write(m_culmative_stream, m_culmative_stream_len * sizeof(char));
	if (m_image_stream_len > 0)
		outF.write(m_image_stream, m_image_stream_len * sizeof(char));
	outF.close();
	m_num_caches_written++;
	return true;
}
bool PyTrans::cacheHeader() {
	if(m_cache_header_len<1)
		return false;
	return sendCacheHeader(
		m_cache_header,
		m_num_grid_headers_cached,
		m_num_image_headers_cached,
		n_PyTrans::data_type_double_code,
		(m_culmative_stream_len + m_image_stream_len)
	);
}
bool PyTrans::sendGrid(
	char ch_stream[],
	const double grid_vals[], 
	int grid_width, 
	int grid_height, 
	int32_t label_code, 
	int32_t axis_code, 
	int32_t start_end_code, 
	int number_of_loops, 
	int jacobi_frame,
	int ch_stream_offset,
	int ch_stream_len
) {
	if (ch_stream == nullptr)
		return false;
	if (ch_stream_len >= 0)
		if (!(getGridStreamLen(grid_width, grid_height) == ch_stream_len))
			return false;
	sendHeader(
		ch_stream, 
		n_PyTrans::grid_code, 
		label_code, 
		axis_code, 
		start_end_code, 
		number_of_loops, 
		jacobi_frame, 
		grid_width, 
		grid_height, 
		1,
		0.f,
		0.f,
		ch_stream_offset);
	int len = grid_width * grid_height;
	int data_stream_offset = m_header_len + ch_stream_offset;
	return n_PyTrans::streamDtoC(len, grid_vals, ch_stream, data_stream_offset);
}
int PyTrans::getGridStreamLen(int grid_width, int grid_height) {
	int d_len = grid_width * grid_height;
	int char_len = n_PyTrans::sizeCharFromL(d_len);
	return char_len + m_header_len;
}
bool PyTrans::sendImage(
	char ch_stream[],
	const unsigned char img_vals[],
	int img_width,
	int img_height,
	int32_t label_code,
	int32_t axis_code,
	int32_t start_end_code,
	int number_of_loops,
	float max,
	float min,
	int expansion,
	int ch_stream_offset,
	int ch_stream_len
) {
	if (ch_stream == nullptr)
		return false;
	int image_size = img_width * img_height * 3;
	if (ch_stream_len >= 0)
		if (!((image_size + m_header_len) == ch_stream_len))
			return false;
	bool retval= sendHeader(
		ch_stream, 
		n_PyTrans::img_code, 
		label_code, 
		axis_code, 
		start_end_code, 
		number_of_loops, 
		0, 
		img_width, 
		img_height, 
		expansion, 
		max, 
		min, 
		ch_stream_offset);
	if (retval) {
		int offset = ch_stream_offset + m_header_len;
		for(int i=0; i<image_size; i++){
			ch_stream[offset + i] = static_cast<char>(img_vals[i]);
		}
	}
	return retval;
}
bool PyTrans::sendCacheHeader(
	char ch_stream[],
	int number_of_stream_headers,
	int number_of_image_headers,
	int data_type_code,
	int cache_size_in_bytes,
	int ch_stream_offset,
	int ch_stream_len
) {
	if (ch_stream_len >= 0 && ch_stream_len > m_cache_header_len)
		return false;
	int32_t cache_header_I[n_PyTrans::cache_header_len];
	cache_header_I[0] = m_cache_header_len;
	cache_header_I[1] = m_header_len;
	cache_header_I[2] = number_of_stream_headers+number_of_image_headers;
	cache_header_I[3] = number_of_stream_headers;
	cache_header_I[4] = number_of_image_headers;
	cache_header_I[5] = data_type_code;
	cache_header_I[6] = cache_size_in_bytes;
	return n_PyTrans::streamItoC(n_PyTrans::cache_header_len, cache_header_I, ch_stream, ch_stream_offset);
}
bool PyTrans::sendHeader(
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
	int ch_stream_offset,
	int ch_stream_len
)
{
	if (ch_stream_len >= 0 && ch_stream_len < m_header_len)
		return false;
	int32_t header_I_stream[PYTRANS_NUMHEADER_FIELDS];
	header_I_stream[0] = type_code;
	header_I_stream[1] = label_code;
	header_I_stream[2] = axis_code;
	header_I_stream[3] = start_end_code;
	header_I_stream[4] = (int32_t)number_of_loops;
	header_I_stream[5] = (int32_t)jacobi_frame;
	header_I_stream[6] = (int32_t)width;
	header_I_stream[7] = (int32_t)height;
	header_I_stream[8] = (int32_t)expansion;
	header_I_stream[9] = n_PyTrans::Fti(max);
	header_I_stream[10] = n_PyTrans::Fti(min);
	return n_PyTrans::streamItoC(PYTRANS_NUMHEADER_FIELDS, header_I_stream, ch_stream, ch_stream_offset);
}