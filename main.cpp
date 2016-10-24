#include <iostream>
#include "MovingAvgFilter.h"
#include "DSPFunctions.hpp"
#include "MappingFunctions.hpp"
#include <Eigen>

using namespace Eigen;
using namespace std;

typedef std::unique_ptr<Ipp32f, void (*)(void *)> unique_Ipp32f; // typedef unique_ptr for data use

// Tester for the DSP Functions
int main()
{
    //-------------------------- MovingAvgFilter Demo -------------------------------------------//
    // Typical matlab use case : dataout = filter(ones(1, smoothWnSize)/smoothWnSize, 1, datain);
    // Where : smoothWnSize = 11
    cout << "[Moving Average Filter Demo]" << endl;
    int smoothWnSize = 11;
    SMA movAvgFilter(smoothWnSize);
    int Naz = 16;
    unique_Ipp32f dataIn(ippsMalloc_32f(Naz), ippsFree);
    for (int i = 0; i < Naz; i++)
        dataIn.get()[i] = 1.0f + i * 0.2f;

    unique_Ipp32f dataOut(ippsMalloc_32f(Naz), ippsFree);
    movAvgFilter.movAvg(dataIn.get(), Naz, 1.0f / smoothWnSize,
                        dataOut.get());

    for (int i = 0; i < Naz; i++)
    {
        std::cout << dataOut.get()[i] << " " << std::endl;
    }

    // ------------------ findnfftinteger demo -> uses the function primeFactors ---------------//
    cout << "[Find FFT Integers Demo]" << endl;
    std::cout << "Finding fftintegers..." << std::endl;
    std::cout << findNFFT(12755) << std::endl;
    std::cout << findNFFT(18000) << std::endl;
    std::cout << findNFFT(13768) << std::endl;
    std::cout << findNFFT(3772) << std::endl;
    std::cout << findNFFT(4000) << std::endl;

    //-------------------------------- Diff demo -----------------------------------------------//
    cout << "[Diff Demo]" << endl;
    int noOfElements = 5;
    unique_Ipp32f diffIn(ippsMalloc_32f(noOfElements), ippsFree);
    unique_Ipp32f diffOut(ippsMalloc_32f(noOfElements - 1), ippsFree);

    for (int i = 0; i < noOfElements; i++)
    {
        diffIn.get()[i] = i;
        cout << diffIn.get()[i] << " ";
    }
    cout << endl;

    diff(diffIn.get(), diffOut.get(), 5);

    for (int i = 0; i < noOfElements - 1; i++)
    {
        cout << diffOut.get()[i] << " ";
    }
    cout << endl;
}
