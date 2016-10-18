### DSP Functions ###

Here are some C++ versions of Matlab DSP Functions, so that you can be saved the trouble of coding them from scratch. I will add more as time goes by.

### Required : ###
- Intel Performance Primitives 9.0 
- C++11 and above support in your compiler

### Nice to Have : ###
- CMake (all compilation parameters are specified in the CMake, so you may modify to suit your environment)

### List of Functions: ###

*Matlab command :*
- filter

*C++ :*
- MovingAvgFilter

<br>
*Matlab command :*
- circshift

*C++ :*
- circshift1D_OP 
- circshift1D_IP 
- circshift (2D)

<br>
*Matlab command :*
- ifftshift 
- fftshift

*C++ :*
- ifftshift1D 
- fftshift1D 
- ifftshift2D 
- fftshift2D

<br>
*Matlab command :*
- sort(in, 'ascend') 
- sort(in, 'descend')

*C++ :*
- sort_indexes

<br>
*Matlab command :*
- cumsum

*C++ :*
- cumulativeSum


<br>
*Matlab command :*
- J:D:K giving a range of elements with D difference between adjacent elements

*C++ :*
- colonRangeVec

<br>
*Matlab command :*
- detrend (linear only)

*C++ :*
- detrend_IP
- detrend_OP

<br>
*Taylor Windowing Function (not matlab in-built but IRL internal function) :*
- taylor_wt

*C++ :*
- Common_TaylorWin

<br>

*Matlab command :*
- factor

*C++ :*
- primeFactor

<br>
*find closest and best fft value for input integer
*C++ :*
- findFFT

<br>
*round to number of significant digits
*C++ :*
- round_to_digits
<br>
*Matlab command :*
- diff

*C++ :*
- diff
<br>
*Matlab command :*
- geodetic2geocentric, geocentric2lvRotationMatrix, geocentric2lv
*C++ :*
- geodetic2geocentric, geocentric2lvRotationMatrix, geocentric2lv


Please see *CMakeList.txt* for the exact build parameters and libraries you need.
Note that geocentric2lv requires Eigen for the matrix computations.
If you have CentOS7,
> yum install eigen3-devel.noarch
> yum install eigen3-doc.noarch

If you don't have Eigen, and want to compile, just remove all references
to the functions in MappingFunctions and from the CMakeLists.txt file.

If this has helped you, kindly help to star this project as a token appreciation. Thanks!