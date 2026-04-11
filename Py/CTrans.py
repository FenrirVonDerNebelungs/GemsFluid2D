from pkgutil import get_data
import struct
from threading import setprofile_all_threads
import numpy as np

#header defs
#1st 32 bits
header1 = {
     0x01 : 'grid',
     0x02 : 'img' 
}
#2nd 32 bits
header2 = {
    0x01: 'U',
    0x02: 'Force',
    0x03: 'relPos',
    0x04: 'W',
    0x05: 'DivW',
    0x06: 'P',
    0x07: 'gradP',
    0x08: 'jacobi_frame',
    0x09: 'DivU'
}
#3rd 32 bits
header3 = {
    0x01: 'X',
    0x02: 'Y',
    0x03: 'Scalar',
    0x04: 'X_exp',
    0x05: 'Y_exp'
    }
#4th 32 bits
header4 = {
    0x00: 'mid_frame',
    0x01 : 'start_frame',
    0x02: 'after_advection',
    0x03: 'after_force',
    0x04: 'end_frame'
    }

cache_header = {
    0x01: 'float',
    0x02: 'double'
}

def getGridWidth(header):
    width = header[6]
    return int(width)

def getGridHeight(header):
    height = header[7]
    return int(height)

def getExpFactor(header):
    exp_factor = header[8]
    return int(exp_factor)

def getDataLen(header): # gets data length in bytes from header, assuming doubles
    if np.size(header)<1:
        return 0
    width = getGridWidth(header)
    height = getGridHeight(header)
    exp_factor = getExpFactor(header)
    d_len = width*height*exp_factor
    return 8*d_len
 
def getDataLabel(header):
    if not header:
        return 0
    return header[1]

def getDataAxis(header):
    return header[2]

def getDataStartEndCode(header):
    return header[3]

def getJacobiFrame(header):
    return header[5]

def getExpansionFactor( header):
    return header[8]

#file header
def getCacheHeaderLenFromFileHeader(header): #returns cache header len in bytes
    return int(header[1])

def getNumberOfCaches(header):
    return int(header[3])

#cache header
def getStreamHeaderLen(header): #getst stream header len in bytes from cache header
    return int(header[1])

def getNumberHeadersPerCache(header):
    return int(header[2])

def getCacheSizeInBytes(header): #size of cache after header
    return int(header[6])

def getDataTypeCode(header):
    raw_type_code = int(header[5])
    return cache_header[raw_type_code]


def Cti(ch4):
    I = int.from_bytes(ch4, "little")
    return I

def Ctl(ch8):
    Long_I = int.from_bytes(ch8, "little")
    return Long_I

def Ctd(ch8):
    double_val = struct.unpack('<d',ch8) #unpack with little-endian ch[7] is the largest value
    return double_val[0]

def Ctf(ch4):
    float_val =struct.unpack('<f',ch4) 
    return float_val

def streamCtoI(len_stream_in, byte_stream_in, stream_in_offset=0):
    ch4 = bytearray(4)
    stream_out_len = int(len_stream_in/4)
    if stream_out_len<1:
        return np.empty(shape=(0,), dtype=np.int32)
    I_stream_out = np.zeros(stream_out_len, dtype=np.int32)
    stream_out_i=0
    for i in range (0, len_stream_in, 4):
        stream_i = stream_in_offset+i
        for ch4_i in range(0,4):
            ch4[ch4_i]=byte_stream_in[stream_i+ch4_i]
        I = Cti(ch4)
        I_stream_out[stream_out_i]=I
        stream_out_i+=1
    return I_stream_out

def streamCtoL(len_stream_in, byte_stream_in, stream_in_offset=0):  #chars are bytes
    ch8 = bytearray(8)
    stream_out_len = len_stream_in/8
    if (stream_out_len<1):
        return np.empty(shape=(0,), dtype=np.int64)
    L_stream_out = np.zeros(stream_out_len, dtype=np.int64)
    stream_out_i=0
    for i in range(0, len_stream_in, 8):
        stream_i=stream_in_offset+i
        for ch8_i in range(0,8):
            ch8[ch8_i]=byte_stream_in[stream_i+ch8_i]
            L_out = Ctl(ch8)
            L_stream_out[stream_out_i]=L_out
            stream_out_i+=1
    return L_stream_out

def streamCtoD(len_stream_in, byte_stream_in, stream_in_offset):
    ch8=bytearray(8)
    stream_out_len = int(len_stream_in/8)
    if stream_out_len<1:
        return np.empty(shape=(0,), dtype=np.float64)
    D_stream_out = np.zeros(stream_out_len, dtype=np.float64)
    stream_out_i = 0
    for i in range(0, len_stream_in, 8):
        stream_i = stream_in_offset+i
        for ch8_i in range(0,8):
            ch8[ch8_i] = byte_stream_in[stream_i+ch8_i]
        D_out = Ctd(ch8)
        D_stream_out[stream_out_i] = D_out
        stream_out_i+=1
    return D_stream_out

def streamCtoF(len_stream_in, byte_stream_in, stream_in_offset):
    ch4=bytearray(4)
    stream_out_len = int(len_stream_in/4)
    if stream_out_len<1:
        return np.empty(shape=(0,), dtype=np.float32)
    stream_out = np.zeros(stream_out_len, dtype=np.float32)
    stream_out_i = 0
    for i in range(0, len_stream_in, 4):
        stream_i = stream_in_offset+i
        for ch4_i in range(0,4):
            ch4[ch4_i] = byte_stream_in[stream_i+ch4_i]
        F_out = Ctf(ch4)
        stream_out[stream_out_i] = F_out
        stream_out_i+=1
    return stream_out

class CTrans:
    
    def __init__(self, filename): #frameSize in bytes
        self.filename = filename
        self.file_header_len=0
        file_header_stream = self.readFileHeader(filename) #sets file_header_len
        self.file_offset = self.file_header_len
        file_header_num_fields = self.file_header_len/4
        file_header = streamCtoI(file_header_num_fields,file_header_stream)
        self.cache_header_len= getCacheHeaderLenFromFileHeader(file_header) #len in byte
        self.cache_header_len_I = self.cache_header_len/4 #len in int32
        self.num_cache_headers = getNumberOfCaches(file_header) #number of caches written to the file
        self.num_caches_read=0
        self.cache_data_type = ''
        self.header_len=0 #stream header len
        self.num_frame_headers = 0
        self.frameSize = 0
        self.frameStream = []

    def readBinaryFile(self, start_offset): #read cache
        retVal=True
        try:
            with open(self.filename, 'rb') as f_in:
                f_in.seek(start_offset)
                cache_header_stream = f_in.read(self.cache_header_len)
                cache_header = streamCtoI(self.cache_header_len_I,cache_header_stream)
                self.header_len = getStreamHeaderLen(cache_header)
                self.num_frame_headers = getNumberHeadersPerCache(cache_header)
                self.frameSize = getCacheSizeInBytes(cache_header)
                
                self.frameStream = f_in.read(self.frameSize)
                self.cache_data_type = getDataTypeCode(cache_header)
                self.num_caches_read +=1
        except FileNotFoundError:
            print(f"Error: could not find file")
            retVal=False
        except Exception as e:
            print(f"Error while reading file: {e}")
            retVal=False
        return retVal

    def readFileHeader(self, filename):
        header_stream = []
        try:
            with open(filename, 'rb') as f_in:
                header_stream_1st = f_in.read(1)
                len_of_header = streamCtoI(1,header_stream_1st)
                self.file_header_len=len_of_header
                header_stream = f_in.read(self.file_header_len-1)
        except FileNotFoundError:
            print(f"Error: could not find file")
            retVal=False
        except Exception as e:
            print(f"Error while reading file: {e}")
            retVal=False
        return header_stream
    
    def readHeader(self, stream_in, stream_offset):
        header_Is = streamCtoI(self.header_len, stream_in, stream_offset)
        return header_Is

    def readGridStream(self, stream_in, stream_offset):
        header = self.readHeader(stream_in, stream_offset)
        data_len = getDataLen(header)
        data = []
        if self.cache_data_type=='double':
            data = streamCtoD(data_len, stream_in, stream_offset+self.header_len)
        elif self.cache_data_type == 'float':
            data = streamCtoF(data_len, stream_in, stream_offset+self.header_len)
        return header, data

    def readDStream(self, stream_in, stream_len, stream_offset):
        char_stream_len = 8*stream_len
        data = streamCtoD(char_stream_len, stream_in, stream_offset)
        return data

    def readFrameStream(self):
        stream_offset=0
        header_stack=[]
        data_stack=[]
        for header_cnt in range(self.num_frame_headers):
            header, data = self.readGridStream(self.frameStream, stream_offset)
            header_stack.append(header)
            data_stack.append(data)
            data_len = getDataLen(header)+self.header_len
            stream_offset=stream_offset+data_len
        return header_stack, data_stack, stream_offset

    def readFrame(self):
        if self.num_caches_read>=self.num_cache_headers:
            return [], [], 0
        if not self.readBinaryFile(self.file_offset):
            return [], [], 0
        header_stack, data_stack, stream_offset = self.readFrameStream()
        new_file_offset = self.file_offset + self.cache_header_len + self.frameSize
        self.file_offset = new_file_offset
        return header_stack, data_stack, new_file_offset

    def readDdat(self, len_in, file_offset):
        if not self.readBinaryFile(file_offset):
            return [], 0
        dat = self.readDStream(self.frameStream, len_in, file_offset)
        new_file_offset=file_offset+len_in
        return dat, new_file_offset