/*
 * DSPFunctions.hpp
 *
 *  Created on: Jan 25, 2016
 *      Author: pier
 */
#ifndef DSPFUNCTIONS_HPP_
#define DSPFUNCTIONS_HPP_

#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <memory>
#include "ipp.h"

struct xyz
{
    double x;
    double y;
    double z;
} typedef xyz;

struct lla
{
    double lat;
    double lon;
    double alt;
} typedef lla;

template<typename T>
std::vector<std::size_t> sort_indexes(const std::vector<T> &v, bool ascending);

template<typename T>
std::vector<T> subVector(const std::vector<T> &vecIn, int firstIdx, int lastIdx);

template<typename T>
void printVec(const std::vector<T> &vecIn);

template<typename T>
void printArray(T *array, std::size_t size);

template<typename T>
std::size_t indexOfLargestElement(T arr[], std::size_t size, T *largestVal);

template<typename T>
std::size_t indexOfLargestElement(T arr[], std::size_t size);

template<typename T>
void circshift(T *in, T *out, int xdim, int ydim, int xshift, int yshift);

template<typename T>
void circshift1D_OP(T *in, T *out, int ydim, int yshift);

template<typename T>
void ifftshift1D(T *in, T *out, int ydim);

template<typename T>
void fftshift1D(T *in, T *out, int ydim);

template<typename T>
void ifftshift2D(T *in, T *out, int xdim, int ydim);

template<typename T>
void fftshift2D(T *in, T *out, int xdim, int ydim);

template<typename T>
void detrend(T *y, int m);

template<typename T>
void matrixTranspose(T *a, int rows, int cols, T *b);

void subArray(const Ipp32f *arrayIn, int start, int end, Ipp32f *subArrayOut);

void cumulativeSum(Ipp32f *array, int size);

void printArrayIpp32fc(Ipp32fc *array, std::size_t size);

void Common_TaylorWin(float *wt, int len);

void simple_transpose_32fc(Ipp32fc *src, Ipp32fc *dst, int nrows, int ncols);

std::vector<int> primeFactors(int n);

int findNFFT(int n);

/**
 *	Index sort
 *
 *	Returns index of data after data is sorted.
 *	Used in Phase Gradient Algorithm.
 *
 *  Inputs:
 *		v - vector of values to be sorted
 *
 *  Outputs:
 *      Returns vector of indices (std::size_t type) of original data elements, sorted.
 *
 */
template<typename T>
std::vector<std::size_t> sort_indexes(const std::vector<T> &v, bool ascending)
{

    // initialize original index locations
    std::vector<std::size_t> idx(v.size());
    for (std::size_t i = 0; i != idx.size(); ++i)
        idx[i] = i;

    if (ascending)
    {
        sort(idx.begin(), idx.end(), [&v](std::size_t i1, std::size_t i2)
        { return v[i1] < v[i2]; });
    }
    else
    {
        sort(idx.begin(), idx.end(), [&v](std::size_t i1, std::size_t i2)
        { return v[i1] > v[i2]; }); //-- descending order
    }

    return idx;
}

/**
 *	subVector
 *
 *	Returns subset of input vector given by firstIdx and lastIdx.
 *
 *  Inputs:
 *		vecIn - Input vector
 *		firstIdx - location of beginning of subset vector
 *		lastIdx - location of end of subset vector
 *
 *  Outputs:
 *      Returns a subset of the input vector.
 *
 */
template<typename T>
std::vector<T> subVector(const std::vector<T> &vecIn, int firstIdx, int lastIdx)
{
    typename std::vector<T>::const_iterator first = vecIn.begin() + firstIdx; // matlab: rgBin0(1:NrgBin1);
    typename std::vector<T>::const_iterator last = vecIn.begin() + lastIdx;
    typename std::vector<T> tempVec(first, last);
    return tempVec;
}

template<typename T>
void printVec(const std::vector<T> &vecIn)
{
    std::cout << "Vector size : " << vecIn.size() << std::endl;
    std::cout << "[";
    for (unsigned int i = 0; i < vecIn.size(); i++)
        std::cout << vecIn[i] << " ";
    std::cout << "]" << std::endl;
}

template<typename T>
void printArray(T *array, std::size_t size)
{
    std::cout << "Array size : " << size << std::endl;
    std::cout << "[";
    std::cout << std::scientific;
    for (unsigned int i = 0; i < size; i++)
        std::cout << i << ": " << array[i] << " ";
    std::cout << "]" << std::endl;
}

/**
 *	indexOfLargestElement
 *
 *	Returns the index of the largest element in the array.
 *
 *  Inputs:
 *		arr[] - input array
 *		size  - size of array
 *		largestVal - pointer to the address to store the largest value
 *
 *	Outputs:
 *		Returns largestIndexValue
 */
template<typename T>
std::size_t indexOfLargestElement(T arr[], std::size_t size, T *largestVal)
{
    std::size_t largestIndex = 0;
    for (std::size_t index = largestIndex; index < size; index++)
    {
        if (arr[largestIndex] < arr[index])
        {
            largestIndex = index;
        }
    }

    *largestVal = arr[largestIndex];
    return largestIndex;
}

/**
 *	indexOfLargestElement
 *
 *	Returns the index of the largest element in the array.
 *
 *  Inputs:
 *		arr[] - input array
 *		size  - size of array
 *		largestVal - pointer to the address to store the largest value
 *
 *	Outputs:
 *		Returns largestIndexValue
 *
 *	Note: Overloaded for cases where you don't need the actual value, just the index.
 */
template<typename T>
std::size_t indexOfLargestElement(T arr[], std::size_t size)
{
    std::size_t largestIndex = 0;
    for (std::size_t index = largestIndex; index < size; index++)
    {
        if (arr[largestIndex] < arr[index])
        {
            largestIndex = index;
        }
    }

    return largestIndex;
}

/**
 *	circshift
 *
 *	Does a 2D circular shift like the Matlab command.
 *
 *  Inputs:
 *		out* - pointer to a buffer for the result
 *		in*	 - pointer to the input data
 *		xdim - size of the x-dimension
 *		ydim - size of the y-dimension
 *		xshift - shifting to be done along x-dimension
 *		yshift - shifting to be done along y-dimension
 *
 *  Can be further optimized using std::rotate
 */
template<typename T>
inline void circshift(T *in, T *out, int xdim, int ydim, int xshift, int yshift)
{
    if (xshift == 0 && yshift == 0)
    {
        out = in; //-- no change
        return;
    }

    for (int i = 0; i < xdim; i++)
    {
        int ii = (i + xshift) % xdim;
        if (ii < 0)
            ii = xdim + ii;
        for (int j = 0; j < ydim; j++)
        {
            int jj = (j + yshift) % ydim;
            if (jj < 0)
                jj = ydim + jj;
            out[ii * ydim + jj] = in[i * ydim + j];
        }
    }
}

/**
 * Does 1D Circshift (in-place)
 *
 * @param in Input array of values, circshift done directly on this
 * @param ydim Length of the array
 * @param yshift Amount to be shifted (+ve is shift right, -ve is shift left)
 */
template<typename T>
inline void circshift1D_IP(T *in, int ydim, int yshift)
{
    if (yshift == 0)
        return;

    if (yshift > 0) // shift right
    {
        //std::rotate(&in[0], &in[ydim - yshift - 1], &in[ydim - 1]);
        std::rotate(in, in + (ydim - yshift), in + ydim);
    }
    else if (yshift < 0) // shift left
    {
        yshift = abs(yshift);
        //std::rotate(&in[0], &in[yshift], &in[ydim - 1]);
        std::rotate(in, in + yshift, in + ydim);
    }

    return;
}


/**
 * Does 1D Circshift (out-of-place)
 *
 * @param in Input array of values
 * @param out Circshifted array of values
 * @param ydim Length of the array
 * @param yshift Amount to be shifted (+ve is shift right, -ve is shift left)
 */
template<typename T>
inline void circshift1D_OP(T *in, T *out, int ydim, int yshift)
{
    if (yshift == 0)
    {
        out = in; //-- no change
        return;
    }

    memcpy(out, in, ydim * sizeof(T));

    if (yshift > 0) // shift right
    {
        //std::rotate(&out[0], &out[ydim - yshift], &out[ydim]); // TODO check indices may be ydim-yshift
        std::rotate(out, out + (ydim - yshift), out + ydim); // C++ idiom: out + ydim is not used, out + ydim -1 is referenced
    }
    else if (yshift < 0) // shift left
    {
        yshift = abs(yshift);
        //std::rotate(&out[0], &out[yshift], &out[ydim - 1]);
        std::rotate(out, out + yshift, out + ydim); // TODO check
    }

    return;

//    for (int j = 0; j < ydim; j++)
//    {
//        int jj = (j + yshift) % ydim;
//        if (jj < 0)
//            jj = ydim + jj;
//        out[jj] = in[j];
//    }
}


/**
 * Does 1D ifftshift
 * Note: T* out must already by memory allocated!!
 */
template<typename T>
inline void ifftshift1D(T *in, T *out, int ydim)
{
    //-- (ydim & 1)==0
    int pivot = (ydim % 2 == 0) ? (ydim / 2) : ((ydim + 1) / 2);
    //circshift1D(in, out, ydim, shiftBy);

    int rightHalf = ydim-pivot;
    int leftHalf = pivot;
    memcpy(out, in+(pivot), sizeof(T)*rightHalf);
    memcpy(out+rightHalf, in, sizeof(T)*leftHalf);
}

/**
 * Does 1D fftshift
 * Note: T* out must already by memory allocated!!
 */
template<typename T>
inline void fftshift1D(T *in, T *out, int ydim)
{
    int pivot = (ydim % 2 == 0) ? (ydim / 2) : ((ydim - 1) / 2);
    //circshift1D(in, out, ydim, shiftBy);
    int rightHalf = ydim-pivot;
    int leftHalf = pivot;
    memcpy(out, in+(pivot), sizeof(T)*rightHalf);
    memcpy(out+rightHalf, in, sizeof(T)*leftHalf);
}

/**
 * Slow due to the circshift, but works.
 */
template<typename T>
inline void ifftshift2D(T *in, T *out, int xdim, int ydim)
{
    int shiftYBy = (ydim % 2 == 0) ? (ydim / 2) : ((ydim + 1) / 2);
    int shiftXBy = (xdim % 2 == 0) ? (xdim / 2) : ((xdim + 1) / 2);
    circshift(in, out, xdim, ydim, shiftXBy, shiftYBy);
}

/**
 * Slow due to the circshift, but works
 */
template<typename T>
inline void fftshift2D(T *in, T *out, int xdim, int ydim)
{
    int shiftYBy = (ydim % 2 == 0) ? (ydim / 2) : ((ydim - 1) / 2);
    int shiftXBy = (xdim % 2 == 0) ? (xdim / 2) : ((xdim - 1) / 2);
    circshift(in, out, xdim, ydim, shiftXBy, shiftYBy);
}

/************************************************************************************
Function    : void detrend_IP(T *y, T *x, int m)
Description : Remove the linear trend of the input floating point data. Note that this
              will initialize a work buffer inside the function. So if you are calling
              this many, many times, create your work buffer in the calling scope and call
              detrend(T *y, T*x, int m) instead to avoid initializing memory over and over
              again.
Inputs      : y - Floating point input data
              m - Input data length
Outputs     : y - Data with linear trend removed
Copyright   : DSO National Laboratories
History     : 01/02/2008, TCK, Adapted from HYC code
              01/12/2008, TCK, Added in return value
              25/01/2016, Pier, Changed into template type, removed need for work buffer
*************************************************************************************/
template<typename T>
void detrend_IP(T *y, int m)
{
    T xmean, ymean;
    int i;
    T temp;
    T Sxy;
    T Sxx;

    T grad;
    T yint;

    std::unique_ptr<T[]> x(new T[m]);

    /********************************
    Set the X axis Liner Values
    *********************************/
    for (i = 0; i < m; i++)
        x[i] = i;

    /********************************
    Calculate the mean of x and y
    *********************************/
    xmean = 0;
    ymean = 0;
    for (i = 0; i < m; i++)
    {
        xmean += x[i];
        ymean += y[i];
    }
    xmean /= m;
    ymean /= m;

    /********************************
    Calculate Covariance
    *********************************/
    temp = 0;
    for (i = 0; i < m; i++)
        temp += x[i] * y[i];
    Sxy = temp / m - xmean * ymean;

    temp = 0;
    for (i = 0; i < m; i++)
        temp += x[i] * x[i];
    Sxx = temp / m - xmean * xmean;

    /********************************
    Calculate Gradient and Y intercept
    *********************************/
    grad = Sxy / Sxx;
    yint = -grad * xmean + ymean;

    /********************************
    Removing Linear Trend
    *********************************/
    for (i = 0; i < m; i++)
        y[i] = y[i] - (grad * i + yint);

}


/************************************************************************************
Function    : void detrend_OP(T *y, T *x, int m)
Description : Remove the linear trend of the input floating point data
Inputs      : y - Floating point input data
              x - Work buffer (must be initialized in calling scope!)
              m - Input data length
Outputs     : y - Data with linear trend removed
Copyright   : DSO National Laboratories
History     : 01/02/2008, TCK, Adapted from HYC code
              01/12/2008, TCK, Added in return value
              25/01/2016, Pier, Changed into template type
*************************************************************************************/
template<typename T>
void detrend_OP(T *y, T*x, int m)
{
    T xmean, ymean;
    int i;
    T temp;
    T Sxy;
    T Sxx;

    T grad;
    T yint;

    /********************************
    Set the X axis Liner Values
    *********************************/
    for (i = 0; i < m; i++)
        x[i] = i;

    /********************************
    Calculate the mean of x and y
    *********************************/
    xmean = 0;
    ymean = 0;
    for (i = 0; i < m; i++)
    {
        xmean += x[i];
        ymean += y[i];
    }
    xmean /= m;
    ymean /= m;

    /********************************
    Calculate Covariance
    *********************************/
    temp = 0;
    for (i = 0; i < m; i++)
        temp += x[i] * y[i];
    Sxy = temp / m - xmean * ymean;

    temp = 0;
    for (i = 0; i < m; i++)
        temp += x[i] * x[i];
    Sxx = temp / m - xmean * xmean;

    /********************************
    Calculate Gradient and Y intercept
    *********************************/
    grad = Sxy / Sxx;
    yint = -grad * xmean + ymean;

    /********************************
    Removing Linear Trend
    *********************************/
    for (i = 0; i < m; i++)
        y[i] = y[i] - (grad * i + yint);

}




/**
 * Works but stupid and slow. Take a look at simple_transpose_32fc
 */
template<typename T>
void matrixTranspose(T *in, int rows, int cols, T *out)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            *(out + (j * rows) + i) = *(in + (i * cols) + j); // a[i][j] == a + i * col + j
        }
    }
}

/*!
 * Same as matlab's diff in 1 dimension
 * f X is a vector, then diff(X) returns a vector, one element shorter than X, of differences between adjacent elements:
 *[X(2)-X(1) X(3)-X(2) ... X(n)-X(n-1)]
 *
 * @param in Input vector
 * @param out Output vector
 * @param noOfElements Self-explanatory
 */
template<typename T>
void diff(T *in, T *out, int noOfElements)
{
    for (int i = 0; i < noOfElements - 1; i++)
        out[i] = in[i + 1] - in[i];
}



#endif /* DSPFUNCTIONS_HPP_ */
