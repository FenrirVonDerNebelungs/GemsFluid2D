import numpy as np
from scipy.signal import convolve2d

def getIndex(i,j,width):
    return j*width+i;
def getValueAtIndex(inVal, i, j, width, height):
    if i<0 or j<0:
        return 0.0
    elif i>=width or j>=height:
        return 0.0
    cur_index = getIndex(i,j, width)
    return inVal[cur_index]

def Lap(inVal, width, height, delta_x):
    outVal = [0]*len(inVal)
    for j in range(height):
        for i in range(width):
                el_L = getValueAtIndex(inVal, i-1, j,width,height)
                el_R = getValueAtIndex(inVal, i+1, j,width,height)
                el_T = getValueAtIndex(inVal, i, j+1, width, height)
                el_B = getValueAtIndex(inVal, i, j-1, width, height)
                el_center =getValueAtIndex(inVal,i,j,width,height)
                lap_val = el_L+el_R+el_T+el_B-4*el_center
                lap_val /= (delta_x*delta_x)
                outVal[getIndex(i,j,width)]=lap_val
    return outVal

def filterIndexLen(index):
    return abs(index-2+0.5)

def genGaussFilter(sigma=1.0):
    filter_array = np.zeros((4,4))
    filter_weight=0.0
    for j in range(4):
        y_len = filterIndexLen(j)
        for i in range(4):
            x_len = filterIndexLen(i)
            r_len_exp2 = x_len*x_len + y_len*y_len
            exp_denom = -1.0/(2.0*sigma**2)
            g_val = np.exp(exp_denom*r_len_exp2)
            filter_array[i][j]=g_val
            filter_weight+=g_val
    filter_array_norm = [el/filter_weight for el in filter_array]
    return filter_array_norm
            

def convolFilter(inVal, filter_kernel, width, height):
    ValMat = np.array(inVal).reshape(height,width)
    ValMat_pad = np.pad(ValMat, ((1,1), (1,1)), mode='constant')
    filtered_with_extra = convolve2d(ValMat_pad, filter_kernel, mode='valid')
    filtered_Val = filtered_with_extra[::2,::2]
    return np.ravel(filtered_Val)

def GaussFilter(inVal,width, height):
    filter_array = genGaussFilter()
    return convolFilter(inVal, filter_array, width, height)

def genReducedStack(baseVal, baseWidth, baseHeight, stack_height): #stack height include base stack is generated in reverse
    stack=[]
    stack.append(baseVal)
    max_stack_index = stack_height-1
    filter_kernel = genGaussFilter()
    width = int(baseWidth)
    height = int(baseHeight)
    for i in range(0, max_stack_index):
        filtered_vals = convolFilter(stack[i], filter_kernel, int(width), int(height))
        width/=2
        height/=2
        stack.append(filtered_vals)
    return stack