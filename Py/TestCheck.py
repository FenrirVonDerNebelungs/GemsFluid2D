import math
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
def Div(inVal_X, inVal_Y, width, height, delta_x):
    outVal=[0]*len(inVal_X)
    inv_d = 1.0/(2.0*delta_x)
    for j in range(height):
        for i in range(width):
            el_L = getValueAtIndex(inVal_X, i-1, j,width,height)
            el_R = getValueAtIndex(inVal_X, i+1, j,width,height)
            el_T = getValueAtIndex(inVal_Y, i, j+1, width, height)
            el_B = getValueAtIndex(inVal_Y, i, j-1, width, height)
            dx = el_R-el_L
            dy = el_T-el_B
            dx *=inv_d
            dy *=inv_d
            div_val = dx+dy
            if i>0 and i<(width-1) and j>0 and j<(height-1):
                outVal[getIndex(i,j,width)]=div_val
    return outVal    

def Dx_wide(inVal, width, height, delta_x):
    outVal=[0]*len(inVal)
    wide_inv_d = 1.0/(2.0*delta_x) 
    short_inv_d = 1.0/delta_x
    inv_d=wide_inv_d               
    for j in range(height):
        for i in range(width):
            inv_d=wide_inv_d
            el_L = getValueAtIndex(inVal, i-1, j,width,height)
            el_R = getValueAtIndex(inVal, i+1, j,width,height)
            if i==0:
                el_L = getValueAtIndex(inVal, i,j,width,height)
                inv_d=short_inv_d
            elif i==(width-1):
                el_R = getValueAtIndex(inVal, i,j,width,height)
                inv_d=short_inv_d
            dx = el_R-el_L
            dx *=inv_d
            outVal[getIndex(i,j,width)]=dx
    return outVal

def Dy_wide(inVal, width, height, delta_x):
    outVal=[0]*len(inVal)
    wide_inv_d = 1.0/(2.0*delta_x) 
    short_inv_d = 1.0/delta_x
    inv_d=wide_inv_d               
    for j in range(height):
        for i in range(width):
            inv_d=wide_inv_d
            el_B = getValueAtIndex(inVal, i, j-1,width,height)
            el_T = getValueAtIndex(inVal, i, j+1,width,height)
            if j==0:
                el_B = getValueAtIndex(inVal, i, j,width,height)
                inv_d=short_inv_d
            elif j==(height-1):
                el_T = getValueAtIndex(inVal, i,j,width, height)
                inv_d=short_inv_d
            dy = el_T-el_B
            dy *=inv_d
            outVal[getIndex(i,j,width)]=dy

    return outVal    
                           
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

def subFields(f1, f2):
    outVal = [val1 - val2 for val1, val2 in zip(f1, f2)]
    return outVal

def addFields(f1,f2):
    outVal = [val1 + val2 for val1, val2 in zip(f1,f2)]
    return outVal

def viscous(inVal, width, height, delta_x, delta_t, nu):
    dW_ = Lap(inVal, width, height, delta_x)
    outVal = [x*nu*delta_t for x in dW_]
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

def findij(index, width):
    j = int( math.floor(index/width))
    i = int(index) - j*int(width)
    return i, j

def dispFromIndexDist(i, j, iL, iR, jB, jT, d1_to_L, d2_to_R, d3_to_B, d4_to_T):
    delta_x = 1.0e-3
    x_dist = 0.0
    y_dist = 0.0
    if iL>=i:
        x_dist = float(iL-i) 
        x_dist += d1_to_L 
    elif iR<=i:
        x_dist = float(iR-i) 
        x_dist -= d2_to_R
    
    if jB >=j:
        y_dist = float(jB-j)
        y_dist += d3_to_B
    elif jT<=j:
        y_dist=float(jT-j)
        y_dist-=d4_to_T
    x_dist*=delta_x
    y_dist*=delta_x
    return x_dist, y_dist
        
    
def dispLocFromIndexes(i1_BL, i2_TL, i3_BR, i4_TR, d1_to_L, d2_to_R, d3_to_B, d4_to_T, width, height):
    displacements_x = np.zeros(len(i1_BL))
    displacements_y = np.zeros(len(i1_BL))
    for j in range(height):
        for i in range(width):
            index_ = j*width + i
            if(i1_BL[index_]>0 and i2_TL[index_]>0 and i3_BR[index_]>0 and i4_TR[index_]>0):
                iL, jB = findij(i1_BL[index_],width)
                iR, jT = findij(i4_TR[index_], width) 
                dL = d1_to_L[index_]
                dR = d2_to_R[index_]
                dB = d3_to_B[index_]
                dT = d4_to_T[index_]  
                x_dist,y_dist = dispFromIndexDist(i,j,iL,iR,jB,jT,dL,dR,dB,dT)
                displacements_x[index_]=x_dist
                displacements_y[index_]=y_dist
    return displacements_x, displacements_y

def dispLocFromU(Ux, Uy, width, height):
    delta_x = 1.0e-3
    delta_t = 1.0e-3
    displacements_x = np.zeros(len(Ux))
    displacements_y = np.zeros(len(Uy))
    for j in range(height):
        for i in range(width):
            index_ = j*width+i
            velocity_x = Ux[index_]
            velocity_y = Uy[index_]
            dx = -velocity_x*delta_t
            dy = -velocity_y*delta_t
            displacements_x[index_]=dx
            displacements_y[index_]=dy
    return displacements_x, displacements_y