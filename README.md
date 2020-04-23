# IIR_Butterworth_Filter_Cpp
C++ code to calculate the coefficients of the Butterworth filter


This code calculate the coefficients of the Band-pass, Band-stop, Low-pass and High-pass Butterworth filters. The file IIR_Butterworth.cpp can be used to test the code. 

Each filter function will return a 2 rows x N coefficients 2D vector, where Row 1 = Numerator and Row 2 = Denumerator

1) Band-pass: the function is "std::vector<std::vector<double> > lp2bp(double, double, int)". The first two arguments are the two cut-off frequencies and the last argument is the order;

2) Band-stop: the function is "std::vector<std::vector<double> > lp2bs(double, double, int)". The first two arguments are the two cut-off frequencies and the last argument is the order;

3) Low-pass: the function is "std::vector<std::vector<double> > lp2lp(double, int)". The first argument is the cut-off frequency and the last argument is the order;

4) High-pass: the function is "std::vector<std::vector<double> > lp2hp(double, int)". The first argument is the cut-off frequency and the last argument is the order;

This code has been written following the Matlab code, so the arguments of each function reflect the arguments that you should pass to the equivalent functions in Matlab. I tried to be consistent with the names of the functions, in case someone wants to compare this code with Matlab code. 

The library Armadillo needs to be downloaded and installed ((http://arma.sourceforge.net/download.html)) along with lapack and blas. I have uploaded the lapack and blas libraries that I have used. Please, note that I had to use #define ARMA_DONT_USE_CXX11 to make armadillo library work with C++/CLR in visual studio 2019. If you don't use C++/CLR template, you should be able to remove this line and based on my communication with the guys who created armadillo, the code should run faster. 
