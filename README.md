# IIR_Butterworth_Filter_Cpp
C++ code to calculate the coefficients of the Butterworth filter


This code calculates the coefficients of the Band-pass, Band-stop, Low-pass and High-pass Butterworth filters. The file IIR_Butterworth.cpp can be used to test the code. It also filters the data, but no zero-phase delay is applied. The name space is: IIR_B_F.

Each filter function will return a 2 rows x N coefficients 2D vector, where Row 1 = Numerator and Row 2 = Denumerator. The method "check_stability_iir" can be used to check the stability of the filter. Please, keep in mind that if the filter is unstable, numerical instability leading to numerical overflow might happen when the order selected is extremely high. If that situation occurs, the program might assign a default value of 10^10 at the denominator and depending on the type of exception even at the numerator.

1) Band-pass: the function is "std::vector<std::vector<double> > lp2bp(double, double, int)". The first two arguments are the two normalized cut-off frequencies (f1/NF, f2/NF), where NF is the Nyquist frequency. This means that the cutoff frequencies must be within the interval of (0,1). The last argument is the order. Please, keep in mind that if you enter order_filt = 2, the order of the filter will be 2 * order_filt = 4;

2) Band-stop: the function is "std::vector<std::vector<double> > lp2bs(double, double, int)". The first two arguments are the two normalized cut-off frequencies (f1/NF, f2/NF), where NF is the Nyquist frequency. This means that the cutoff frequencies must be within the interval of (0,1). The last argument is the order. Please, keep in mind that if you enter order_filt = 2, the order of the filter will be 2 * order_filt = 4;

3) Low-pass: the function is "std::vector<std::vector<double> > lp2lp(double, int)". The first argument is the normalized cut-off frequency (f/NF), where NF is the Nyquist frequency. This means that the cutoff frequency must be within the interval of (0,1). The last argument is the order;

4) High-pass: the function is "std::vector<std::vector<double> > lp2hp(double, int)". The first argument is the normalized cut-off frequency (f/NF), where NF is the Nyquist frequency. This means that the cutoff frequency must be within the interval of (0,1). The last argument is the order;

5) Check the stability of the filter: the method is "bool check_stability_iir(std::vector<std::vector<double> >)". The argument is the 2D array containing the filter coefficients. It returns "true" if the filter is stable, "false" if it is unstable. 

6) Filter the data: the method is "std::vector<double> Filter_Data(std::vector<std::vector<double> > coeff_filt, std::vector<double> pre_filt_signal)". The two arguments are the filter coefficients and the signal to be filtered. It returns the filtered signal.

The library Armadillo needs to be downloaded and installed ((http://arma.sourceforge.net/download.html)) along with lapack and blas. I have uploaded the lapack and blas libraries that I have used. Please, note that with the older version of Armadillo, I had to use #define ARMA_DONT_USE_CXX11 to make armadillo library work with C++/CLR in visual studio 2019. If you use the latest version (armadillo-9.880.1), which I would recommend, because it is supposedly faster than the previous one, as the developers told me, you should replace #define ARMA_DONT_USE_CXX11 with #define ARMA_DONT_USE_CXX11_MUTEX. 

If you are running the code in Linux (Ubuntu 20.04), you need to make the following changes in IIR_Butterworth.cpp:
1) Comment out #include <crtdbg.h>
2) Cooment out  _CrtDumpMemoryLeaks() in lines 64, 125, 187 and 246.

You can also remove the line of code #define ARMA_DONT_USE_CXX11_MUTEX (or #define ARMA_DONT_USE_CXX11) and, based on my communication with the guys who created armadillo, the code should run faster. 

Compile the code in the following way: 
  
  1) g++ -c IIR_Butterworth_H.cpp (to generate IIR_Butterworth_H.o)
  2) g++ -ggdb IIR_Butterworth.cpp IIR_Butterworth_H.o -larmadillo -o <Name_Exe_File>

If you have any question and/or want to report bugs, please e-mail me (Ale) at: pressalex@hotmail.com
