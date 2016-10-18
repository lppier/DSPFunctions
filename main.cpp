#include <iostream>
#include "MovingAvgFilter.h"
#include "DSPFunctions.hpp"
#include "MappingFunctions.hpp"

using namespace std;

typedef std::unique_ptr<Ipp32f, void (*)(void *)> unique_Ipp32f; //-- typedef unique_ptr for data use

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


    //------------------------------- Testing Mapping Conversion Functions -----------------------------------//
    cout << "[Geodetic -> Geocentric -> ENU Conversion Demo]" << endl;
    int numNavPulses = 13197;

    FILE *fpRs = fopen("lat_lon_alt.dat", "rb");
    if (fpRs == NULL)
    {
        printf("%s(%d): Data file not found\n", __FILE__, __LINE__);
        return 0;
    }
    double3 *m_Rs = (double3 *) malloc(
            numNavPulses * sizeof(double3));
    fread(m_Rs, sizeof(double3), numNavPulses, fpRs);
    fclose(fpRs);

    double3 *m_Output = (double3 *) malloc(numNavPulses * sizeof(double3));
    double pi = 3.14159;
    double ref_Lat = m_Rs[0].x / 180.0 * pi;
    double ref_Lon = m_Rs[0].y / 180.0 * pi;
    double ref_Alt = m_Rs[0].z;


    // Convert lat-lon to local ENU
    for (int i = 0; i < numNavPulses; i++)
    {
        m_Output[i] = geodetic2geocentric(m_Rs[i].x / 180.0 * pi, m_Rs[i].y / 180.0 * pi, m_Rs[i].z);
        m_Output[i] = geocentric2lv(m_Output[i].x, m_Output[i].y, m_Output[i].z, ref_Lat,
                                    ref_Lon, ref_Alt);
        cout << "enu [" << m_Output[i].x << ", " << m_Output[i].y << ", " << m_Output[i].z << "]" << endl;
    }

    free(m_Rs);
    free(m_Output);

    return 0;
}