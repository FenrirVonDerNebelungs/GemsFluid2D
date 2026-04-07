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

def getGridWidth(header):
    width = header[6]
    return width

def getGridHeight(self, header):
    height = header[7]
    return height

def getExpFactor(self, header):
    exp_factor = header[8]
    return exp_factor

def getDataLen(self, header): # gets data length in bytes from header, assuming doubles
    if not header:
        return 0
    width = self.getGridWidth(self, header)
    height = self.getGridHeight(self, header)
    exp_factor = self.getExpFactor(self, header)
    d_len = width*height*exp_factor
    return 8*d_len
 
def getDataLabel(self, header):
    if not header:
        return 0
    return header[1]

def getDataAxis(self,header):
    return header[2]

def getDataStartEndCode(self, header):
    return header[3]

def getJacobiFrame(self, header):
    return header[5]

def getExpansionFactor(self, header):
    return header[8]

def Cti(ch4):
    I = int.from_bytes(ch4, "little")
    return I

def Ctl(ch8):
    Long_I = int.from_bytes(ch8, "little")
    return Long_I

def Ctd(ch8):
    double_val = struct.unpack('<d',ch8) #unpack with little-endian ch[7] is the largest value
    return double_val

def streamCtoI(len_stream_in, byte_stream_in, stream_in_offset=0):
    ch4 = bytearray(4)
    stream_out_len = len_stream_in/4
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
            stream_out_I+=1
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
    stream_out_len = len_stream_in/8
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

class CTrans:
    num_header_fields=9
    header_len = 9*4 # 9 header fields * len int32
    
    def __init__(self, filename, num_headers_per_frame, frameSize): #frameSize in bytes
        self.filename = filename
        self.num_frame_headers = num_headers_per_frame
        self.frameSize = frameSize
        self.frameStream = bytearray(frameSize)

    def readBinaryFile(self, start_offset):
        retVal=True
        try:
            with open(self.filename, 'rb') as f_in:
                f_in.seek(start_offset)
                self.frameStream = f_in.read(self.frameSize)
        except FileNotFoundError:
            print(f"Error: could not find file")
            retVal=False
        except Exception as e:
            print(f"Error while reading file: {e}")
            retVal=False
        return retVal

    def readHeader(self, stream_in, stream_offset):
        header_Is = streamCtoI(self.header_len, stream_in, stream_offset)
        return header_Is

    def readGridStream(self, stream_in, stream_offset):
        header = self.readHeader(stream_in, stream_offset)
        data_len = getDataLen(header)
        data = streamCtoD(data_len, stream_in, stream_offset+self.header_len)
        return header, data

    def readDStream(self, stream_in, stream_len, stream_offset):
        data = streamCtoD(stream_len, stream_in, stream_offset)
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

    def readFrame(self, file_offset):
        if not self.readBinaryFile(file_offset):
            return [], [], 0
        header_stack, data_stack, stream_offset = self.readFrameStream()
        new_file_offset = file_offset+self.frameSize
        return header_stack, data_stack, new_file_offset

    def readDdat(self, len_in, file_offset):
        if not self.readBinaryFile(file_offset):
            return [], 0
        dat = self.readDStream(self.frameStream, len_in, file_offset)
        new_file_offset=file_offset+len_in
        return dat, new_file_offset