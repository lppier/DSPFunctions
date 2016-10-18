#include "DSPFunctions.hpp"
//#include <ipp.h>

void printArrayIpp32fc(Ipp32fc *array, std::size_t size)
{
    std::cout << "Array size : " << size << std::endl;
    std::cout << "[";
    std::cout << std::scientific;
    for (unsigned int i = 0; i < size; i++)
        std::cout << i << ": " << array[i].re << "," << array[i].im << " " << std::endl;
    std::cout << "]" << std::endl;
    std::cout << " " << std::endl;
}


/**
 *	subArray
 *
 *	Returns subset of input array given by start and end.
 *
 *  Inputs:
 *		arrayIn - Input array
 *		start - location of beginning of subset vector
 *		end - location of end of subset vector
 *		subArrayOut - output subset array
 *
 *  Notes: The output arary is inclusive of arrayIn[start] and arrayIn[end]
 */
void subArray(const Ipp32f *arrayIn, int start, int end, Ipp32f *subArrayOut)
{
    //Ipp32f mag_sum[end-start + 1];
    int count = 0;

    for (int i = start; i <= end; i++)
    {
        subArrayOut[count] = arrayIn[i];
        count++;
    }
}

/**
 *	mean
 *
 *	Returns mean of input array subset given by start and end.
 *
 *  Inputs:
 *		arrayIn - Input array
 *		start - location of beginning of array to calculate mean
 *		end - location of end of array for mean calculation
 *
 *	Outputs:
 *		mean
 *
 *  Notes:
 *      Total is inclusive of array[start] and array[end] when calculating mean
 */
Ipp32f mean(Ipp32f array[], int start, int end)
{

    double total = 0.0;
    for (int i = start; i <= end; i++)
    {
        total += array[i];
//		std::cout << std::scientific;
//		std::cout << i << " : " << array[i] << std::endl;
    }
//	std::cout << std::scientific;
//	std::cout << "total" << " : " << total << std::endl;
    return (Ipp32f) (total / (end - start + 1));
}

/**
 *	stddev
 *
 *	Returns standard deviation of input array subset given by start and end.
 *
 *  Inputs:
 *		arrayIn - Input array
 *		start - location of beginning of array to calculate standard deviation
 *		end - location of end of array for standard deviation calculation
 *
 *	Outputs:
 *		standard deviation
 *
 *  Notes:
 *      Total is inclusive of array[start] and array[end] when calculating standard deviation
 */
Ipp32f stddev(Ipp32f array[], int start, int end)
{

    long double sum = 0.0;
    Ipp32f meanVar = mean(array, start, end);

    for (int j = start; j <= end; j++)
    {
        //sum += pow((array[j]-meanVar), 2);
        //long double mulVal = (array[j]-meanVar) * (array[j]-meanVar);
        long double mulVal = pow((array[j] - meanVar), 2);
        sum += mulVal;
    }


    return (Ipp32f) sqrt((sum / (end - start))); //--  -1 from (end-start+1) (stddev formula has -1)
}

/**
 *	cumulativeSum
 *
 *	Returns the cumulative sum of the elements in the array
 *
 *  Inputs:
 *		arrayIn - Input array
 *		start - location of beginning of array to calculate standard deviation
 *		end - location of end of array for standard deviation calculation
 *
 *	Outputs:
 *		standard deviation
 *
 *  Notes:
 *      Total is inclusive of array[start] and array[end] when calculating standard deviation
 */
void cumulativeSum(Ipp32f *array, int size)
{

    if (size < 0) return;
    cumulativeSum(array, size - 1);
    array[size + 1] += array[size];
    //std::cout << "size[" << size+1 << "] += [ " << size << "]" << std::endl;
}


/**
 *	colonRangeVec
 *
 *	Returns a vector of range values according to the Matlab J:D:K syntax
 *	J:D:K  is the same as [J, J+D, ..., J+m*D] where m = fix((K-J)/D).
 *
 *  Inputs:
 *		startVal - Equivalent of J
 *		granularity - Equivalent of D
 *		endVal - Equivalent of J+m*D
 *
 *	Outputs:
 *		Returns vector of float values
 */
std::vector<double> colonRangeVec(double J, double D, double K)
{
    int m = static_cast<int>(((K - J) / D));
    std::vector<double> v;

    //-- Return empty vector if hit the cases below
    if (D == 0 || ((D > 0) && (J > K)) || ((D < 0) && (J < K)))
    {
        return v;
    }

    //-- Else create the vector
    for (int i = 0; i <= m; i++)
    {
        double newVal = J + i * D;
        v.insert(v.begin()+i, newVal);
    }

    return v;
}


/**
 *	round_to_digits
 *
 *	Rounds a double value to number of significant digits specified by input.
 *
 *  Inputs:
 *  	value - the input value to be rounded
 *  	digits - number of significant digits
 *
 *	Outputs:
 *		Returns double value rounded to the specified number of significant digits
 */
double round_to_digits(double value, int digits)
{
    if (value == 0.0)
        return 0.0;

    double factor = pow(10.0, digits - ceil(log10(fabs(value))));
    return round(value*factor) / factor;
}

//-- Equivalent to what IRL is using for Taylor Window (Ken Yew)
void Common_TaylorWin(float *wt, int len)
{
    /******************************
        Constant Parameters -
        Please ask Ken Yew if u want to change (quit liao)
     ******************************/
    int n = 5;
    double SLL = -35.0;
    const double pi = 4.0f * atan(1.0f);
    /******************************
        Parameters
     ******************************/
    double A, B;
    double sigma2;
    int m, i, k;
    double tmp[5][5];
    double Fm_num[5] = {0};
    double Fm_den[5] = {0};
    double Fm[5] = {0};

    /******************************
        Compute Constants
     ******************************/
    B = pow(10.0, -SLL / 20.0);
    A = log(B + sqrt(B * B - 1)) / pi;
    sigma2 = (double) n * (double) n / (A * A + ((double) n - 0.5) * ((double) n - 0.5));

    /******************************
        Generate Numerator
     ******************************/
    /*lint -e834 */
    for (m = 1; m < n; m++)
    {
        Fm_num[m] = 1.0;
        for (i = 1; i < n; i++)
        {
            tmp[m][i] = 1 - (double) m * (double) m / sigma2 * (1 / (A * A + ((double) i - 0.5) * ((double) i - 0.5)));
            Fm_num[m] = Fm_num[m] * tmp[m][i];
        }
        Fm_num[m] = pow(-1.0, (double) (m + 1)) * Fm_num[m];
    }

    /******************************
        Generate Denominator
     ******************************/
    for (m = 1; m < n; m++)
    {
        Fm_den[m] = 1.0;
        for (i = 1; i < n; i++)
        {
            tmp[m][i] = 1 - ((double) m * (double) m) / ((double) i * (double) i);
            if (i == m)
            {
                tmp[m][i] = 1.0;
            }
            Fm_den[m] = Fm_den[m] * tmp[m][i];
        }
        Fm_den[m] = 2 * Fm_den[m];
    }

    /******************************
        Generate Weights
     ******************************/
    for (m = 1; m < n; m++)
        Fm[m] = Fm_num[m] / Fm_den[m];

    for (k = 0; k < len; k++)
    {
        wt[k] = 0;
        for (m = 1; m < n; m++)
            wt[k] += (float) (Fm[m] *
                              cos((2 * pi * (double) m * ((double) k - (double) len / 2.0 + 0.5)) / (double) len));
        wt[k] = 1 + 2 * wt[k];
    }
}

//-- Transpose Ipp32fc matrix using ippi. (Even though we are using ipps mostly)
//-- Can do this because 16bits * 4 = 64bits == sizeof(Ipp32fc)
//-- Note: If Ipp32f, should use ippiTranspose_8u_C4R (32 bits) == sizeof(Ipp32f)
void simple_transpose_32fc(Ipp32fc *src, Ipp32fc *dst, int nrows, int ncols)
{
    int src_stride = ncols * sizeof(*src);
    int dst_stride = nrows * sizeof(*dst);
    // Note that IPPI uses [col, row] for Roi
    IppiSize srcRoi = {ncols, nrows};
    ippiTranspose_16u_C4R((Ipp16u *) src, src_stride, (Ipp16u *) dst, dst_stride, srcRoi);
}

/**
 * @brief Returns a vector containing the prime factors of n
 *
 * @param [in] The number to find the prime factors for
 * @return
 */
std::vector<int> primeFactors(int n)
{
    std::vector<int> vec;

    while (n % 2 == 0)
    {
        vec.push_back(2);
        n /= 2;
    }

    for (int i = 3; i <= sqrt(n); i += 2)
    {
        while (n % i == 0)
        {
            vec.push_back(i);
            n /= i;
        }
    }

    if (n > 2)
        vec.push_back(n);

//    std::cout << "Prime factors:" << std::endl;
//    for (int j=0; j < vec.size(); j++)
//    {
//        printf("%d ", vec[j]);
//    }
//    printf("\n");
    return vec;
}

/**
 * @brief Used to find the appropriate fft integer for the input n
 * This uses the "formula" (N + D - 1)/D * D
 * Criteria: Output nfft should be a factor of 2,3,5
 *
 * @param [in] Integer to find nfft for
 */
int findNFFT(int n)
{
    std::vector<int> ansPrimes;
    std::vector<int> firstPrimes;

    int d = 0;

    do
    {
        if (n > 2048) d = 512;
        else if (n > 1024) d = 256;
        else if (n > 128) d = 64;
        else if (n > 32) d = 32;
        else if (n > 8) d = 8;
        else d = 2;

        int fn = (n + d - 1) / d * d;
        firstPrimes = primeFactors(fn);

        for (int i = 0; i < firstPrimes.size(); i++)
        {
            if (firstPrimes[i] == 2 || firstPrimes[i] == 3 || firstPrimes[i] == 5)
            {
                ansPrimes.push_back(firstPrimes[i]);
                firstPrimes.erase(firstPrimes.begin() + i);
                i -= 1;
            }
        }

        int newN = 1;
        if (firstPrimes.size() > 0)
        {
            for (int i = 0; i < firstPrimes.size(); i++)
                newN *= firstPrimes[i];
        }

        n = newN;
        firstPrimes = {};

    } while (n != 1); // if n == 1 means that firstPrimes

    int ans = 1;
    for (int i = 0; i < ansPrimes.size(); i++)
        ans *= ansPrimes[i];

    return ans;
}




