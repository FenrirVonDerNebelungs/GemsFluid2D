#pragma once
#ifndef BASE_H
#define BASE_H
#include <stdio.h>
#include <cmath>
#include <cstring>
struct s_WH {
	int width;
	int height;
};
struct s_frame_index {
    int in;
    int out;
};
struct s_jacobi_vars {
    double alpha;
    double beta;
    double rbeta;
};
struct s_four_corners {
    int i1;
    int i2;
    int j1;
    int j2;
};
struct s_force {
    int i;
    int j;
    double Fx_c;
    double Fy_c;
    double R;/*falloff radius*/
    double inv_Rsqrd;/* 1/R^2 */
	bool active;
};

inline void reverseFrameIndex(int& frame_index) {
    frame_index = (frame_index + 1) % 2;
}
inline s_frame_index getFrameIndex(int frame_index) {
	s_frame_index f_indexes;
    f_indexes.in = frame_index;
    f_indexes.out = (f_indexes.in + 1) % 2;
    return f_indexes;
}
inline void setFrameIndexes(s_frame_index& f_indexes, int frame_index) {
    f_indexes.in = frame_index;
    f_indexes.out = (f_indexes.in + 1) % 2;
}
inline void swapFrameIndexes(s_frame_index& f_indexes) {
    f_indexes.in = f_indexes.out;
    f_indexes.out = (f_indexes.in + 1) % 2;
}
inline void swapFramePointers(double* FramePtr[]) {
    double* swap = FramePtr[0];
    FramePtr[0] = FramePtr[1];
    FramePtr[1] = swap;
}
inline void fixFramePointers(double* FramePtr[], s_frame_index& f_indexes, int original_frame_in_index) {
    if(original_frame_in_index != f_indexes.in) {
        swapFramePointers(FramePtr);
	}
}
inline bool inRange_ij_CPU(int i, int j, int width, int height) {
    return (i >= 0 && i < width && j >= 0 && j < height);
}
#endif
