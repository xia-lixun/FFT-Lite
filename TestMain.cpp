/* 
 * File:   TestMain.cpp
 * Author: lixun
 *
 * Created on January 24, 2015, 4:32 PM
 */




//R2 SP Time Elapsed:520692361 ns
//R4 SP Time Elapsed:295550860 ns
//R2 DP Time Elapsed:1156191757 ns
//R4 DP Time Elapsed:703803318 ns
//FFTWf Time Elapsed:85714582 ns
//FFTW Time Elapsed:149044621 ns
//err R2:5.215124e-04             err R4:2.281617e-04
//err R2:2.455753835930862e-12    err R4:2.405840844699303e-13




#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
#include <cstring>
#include <limits>
#include <fftw3.h>


#include "Config.h"

#include "CSort.h"
#include "CSort4.h"

#include "CCoef.h"
#include "CCoefd.h"
#include "CCoef4.h"
#include "CCoef4d.h"

#include "CTrans.h"
#include "CTrans4.h"

#include "CChirpZ.h"
#include "CChirpZd.h"

#include "CFastFourier.h"
#include "CFastBluestein.h"





#define N_TEST_RUN    (1)
//#define SHOW_VALUES
#define TEST_USE_FFTW (1)



using namespace std;







void test_fft_f(int len) {
    
    float * in = new float[len*2];
    float * out = new float[len*2];
    float * ref = new float[len*2];    
    
    srand(time(NULL));
    for(int i = 0; i < len*2; i++) {
        in[i] = (float) ((double)rand() / (double)RAND_MAX - 0.5);        
    }
    
    unique_ptr<CFastFourier> pFFTL(new CFastFourier);
    pFFTL->fft_f(len, in, out);
    
    
    
    
    // reference output
#ifdef TEST_USE_FFTW    
    fftwf_complex *inFFTW, *outFFTW;
    fftwf_plan pFFTW;
    
    inFFTW = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * len);
    outFFTW = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * len);
    pFFTW = fftwf_plan_dft_1d(len, inFFTW, outFFTW, FFTW_FORWARD, FFTW_MEASURE);    
    memcpy(inFFTW, in, sizeof(float)*len*Cx2Re);         
    fftwf_execute(pFFTW);    
    memcpy(ref, outFFTW, sizeof(float)*len*Cx2Re);

    
    //statistics
    
    double error  = 0.0;
    double error_ = 0.0;
    int idx = 0;
    for(int i = 0; i < len*2; i++) {
        if(fabs(out[i] - ref[i]) > error) {
            error = fabs(out[i] - ref[i]);
            error_ = ref[i];
            idx = i;
        }
    }
    cout.precision(std::numeric_limits<double>::digits10);
    cout << scientific << "FFTL/FFTW Forward Float:" << error/error_ << "@" << idx << endl;    

    
    // clear the site
    
    fftwf_destroy_plan(pFFTW);
    fftwf_free(inFFTW);
    fftwf_free(outFFTW);
#endif    
    
    delete[] in;
    delete[] out;
    delete[] ref;    
}





void test_fft_d(int len) {
    
    double * in = new double[len*2];
    double * out = new double[len*2];
    double * ref = new double[len*2];    
    
    srand(time(NULL));
    for(int i = 0; i < len*2; i++) {
        in[i] = (double)rand() / (double)RAND_MAX - 0.5;        
    }
    
    unique_ptr<CFastFourier> pFFTL(new CFastFourier);
    pFFTL->fft_d(len, in, out);
    
    
    
    
    // reference output
#ifdef TEST_USE_FFTW        
    fftw_complex *inFFTW, *outFFTW;
    fftw_plan pFFTW;
    
    inFFTW = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * len);
    outFFTW = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * len);
    pFFTW = fftw_plan_dft_1d(len, inFFTW, outFFTW, FFTW_FORWARD, FFTW_MEASURE);    
    memcpy(inFFTW, in, sizeof(double)*len*Cx2Re);         
    fftw_execute(pFFTW);    
    memcpy(ref, outFFTW, sizeof(double)*len*Cx2Re);
    
    
    //statistics
    
    double error  = 0.0;
    double error_ = 0.0;
    int idx = 0;
    for(int i = 0; i < len*2; i++) {
        if(fabs(out[i] - ref[i]) > error) {
            error = fabs(out[i] - ref[i]);
            error_ = ref[i];
            idx = i;
        }
    }
    cout.precision(std::numeric_limits<double>::digits10);
    cout << scientific << "FFTL/FFTW Forward Double:" << error/error_ << "@" << idx << endl;    
    
    
    // clear the site
    
    fftw_destroy_plan(pFFTW);
    fftw_free(inFFTW);
    fftw_free(outFFTW);
#endif    
    delete[] in;
    delete[] out;
    delete[] ref;    
}







void test_ifft_f(int len) {

    float * in = new float[len*2];
    float * out = new float[len*2];
    float * ref = new float[len*2];    
    
    srand(time(NULL));
    for(int i = 0; i < len*2; i++) {
        in[i] = (float) ((double)rand() / (double)RAND_MAX - 0.5);        
    }
    
    unique_ptr<CFastFourier> pFFTL(new CFastFourier);
    pFFTL->ifft_f(len, in, out);
    
    
    
    
    // reference output
#ifdef TEST_USE_FFTW        
    fftwf_complex *inFFTW, *outFFTW;
    fftwf_plan pFFTW;
    
    inFFTW = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * len);
    outFFTW = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * len);
    pFFTW = fftwf_plan_dft_1d(len, inFFTW, outFFTW, FFTW_BACKWARD, FFTW_MEASURE);    
    memcpy(inFFTW, in, sizeof(float)*len*Cx2Re);         
    fftwf_execute(pFFTW);    
    memcpy(ref, outFFTW, sizeof(float)*len*Cx2Re);
    
    
    //statistics
    
    double error  = 0.0;
    double error_ = 0.0;
    int idx = 0;
    for(int i = 0; i < len*2; i++) {
        if(fabs(out[i] - ref[i]) > error) {
            error = fabs(out[i] - ref[i]);
            error_ = ref[i];
            idx = i;
        }
    }
    cout.precision(std::numeric_limits<double>::digits10);
    cout << scientific << "FFTL/FFTW Backward Float:" << error/error_ << "@" << idx << endl;    
    
    
    // clear the site
    
    fftwf_destroy_plan(pFFTW);
    fftwf_free(inFFTW);
    fftwf_free(outFFTW);
#endif    
    
    delete[] in;
    delete[] out;
    delete[] ref;    
}





void test_ifft_d(int len) {
    
    double * in = new double[len*2];
    double * out = new double[len*2];
    double * ref = new double[len*2];    
    
    srand(time(NULL));
    for(int i = 0; i < len*2; i++) {
        in[i] = (double)rand() / (double)RAND_MAX - 0.5;        
    }
    
    unique_ptr<CFastFourier> pFFTL(new CFastFourier);
    pFFTL->ifft_d(len, in, out);
    
    
    
    
    // reference output
#ifdef TEST_USE_FFTW        
    fftw_complex *inFFTW, *outFFTW;
    fftw_plan pFFTW;
    
    inFFTW = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * len);
    outFFTW = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * len);
    pFFTW = fftw_plan_dft_1d(len, inFFTW, outFFTW, FFTW_BACKWARD, FFTW_MEASURE);    
    memcpy(inFFTW, in, sizeof(double)*len*Cx2Re);         
    fftw_execute(pFFTW);    
    memcpy(ref, outFFTW, sizeof(double)*len*Cx2Re);
    
    
    //statistics
    
    double error  = 0.0;
    double error_ = 0.0;
    int idx = 0;
    for(int i = 0; i < len*2; i++) {
        if(fabs(out[i] - ref[i]) > error) {
            error = fabs(out[i] - ref[i]);
            error_ = ref[i];
            idx = i;
        }
    }
    cout.precision(std::numeric_limits<double>::digits10);
    cout << scientific << "FFTL/FFTW Backward Double:" << error/error_ << "@" << idx << endl;    
    
    
    // clear the site
    
    fftw_destroy_plan(pFFTW);
    fftw_free(inFFTW);
    fftw_free(outFFTW);
#endif 
    
    delete[] in;
    delete[] out;
    delete[] ref;      
}



//FFTL/FFTW Forward Float:-6.001999866355470e-07@12725215
//FFTL/FFTW Forward Float:-4.577024572442971e-07@242414
//FFTL/FFTW Backward Float:4.178242138324048e-07 @ 19789539
//FFTL/FFTW Backward Float:1.223873715814507e-06 @ 2357043
//FFTL/FFTW Backward Float:7.345452595570766e-07 @ 10275056

//FFTL/FFTW Forward Double:2.803627078783935e-14@50326501
//FFTL/FFTW Forward Double:-3.302118814100065e-15@10379231
//FFTL/FFTW Backward Double:3.255823608056726e-15@16070462
//FFTL/FFTW Backward Double:4.229361071603665e-14@32987122








void test_chirpz_f(int len) {
    
    float * in = new float[len*2];
    float * out = new float[len*2];
    float * ref = new float[len*2];    
    
    srand(time(NULL));
    for(int i = 0; i < len*2; i++) {
        in[i] = (float) ((double)rand() / (double)RAND_MAX - 0.5);        
    }
    
    unique_ptr<CFastBluestein> pChirpZ(new CFastBluestein);
    pChirpZ->chirpz_f(len, in, out);
    
#ifdef SHOW_VALUES   
    std::cout << "FFTL:" << std::endl;
    cout.precision(std::numeric_limits<float>::digits10);        
    for(int k = 0; k < len; k++) {
        std::cout << scientific << out[re(k)] << "    "<< out[im(k)] << std::endl;
    }
    std::cout << "FFTW:" << std::endl;    
#endif 
    
    
    
    
    // reference output
#ifdef TEST_USE_FFTW        
    fftwf_complex *inFFTW, *outFFTW;
    fftwf_plan pFFTW;
    
    inFFTW = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * len);
    outFFTW = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * len);
    pFFTW = fftwf_plan_dft_1d(len, inFFTW, outFFTW, FFTW_FORWARD, FFTW_MEASURE);    
    memcpy(inFFTW, in, sizeof(float)*len*Cx2Re);         
    fftwf_execute(pFFTW);    
    memcpy(ref, outFFTW, sizeof(float)*len*Cx2Re);
#ifdef SHOW_VALUES
    cout.precision(std::numeric_limits<float>::digits10);    
    for(int k = 0; k < len; k++) {
        std::cout << scientific << ref[re(k)] << "    "<< ref[im(k)] << std::endl;
    }    
#endif
    //statistics
    
    double error  = 0.0;
    double error_ = 0.0;
    int idx = 0;
    for(int i = 0; i < len*2; i++) {
        if(fabs(out[i] - ref[i]) > error) {
            error = fabs(out[i] - ref[i]);
            error_ = ref[i];
            idx = i;
        }
    }
    cout.precision(std::numeric_limits<double>::digits10);
    cout << scientific << "FFTL/FFTW Forward Float:" << error/error_ << "@" << idx << endl;    
    
    
    // clear the site
    
    fftwf_destroy_plan(pFFTW);
    fftwf_free(inFFTW);
    fftwf_free(outFFTW);
    
#endif   
 
    delete[] in;
    delete[] out;
    delete[] ref;    
}




void test_ichirpz_f(int len) {
    
    float * in = new float[len*2];
    float * out = new float[len*2];
    float * ref = new float[len*2];    
    
    srand(time(NULL));
    for(int i = 0; i < len*2; i++) {
        in[i] = (float) ((double)rand() / (double)RAND_MAX - 0.5);        
    }
    
    unique_ptr<CFastBluestein> pChirpZ(new CFastBluestein);
    pChirpZ->ichirpz_f(len, in, out);
    
#ifdef SHOW_VALUES   
    std::cout << "FFTL:" << std::endl;
    cout.precision(std::numeric_limits<float>::digits10);        
    for(int k = 0; k < len; k++) {
        std::cout << scientific << out[re(k)] << "    "<< out[im(k)] << std::endl;
    }
    std::cout << "FFTW:" << std::endl;    
#endif 
    
    
    
    
    // reference output
#ifdef TEST_USE_FFTW        
    fftwf_complex *inFFTW, *outFFTW;
    fftwf_plan pFFTW;
    
    inFFTW = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * len);
    outFFTW = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * len);
    pFFTW = fftwf_plan_dft_1d(len, inFFTW, outFFTW, FFTW_BACKWARD, FFTW_MEASURE);    
    memcpy(inFFTW, in, sizeof(float)*len*Cx2Re);         
    fftwf_execute(pFFTW);    
    memcpy(ref, outFFTW, sizeof(float)*len*Cx2Re);
#ifdef SHOW_VALUES
    cout.precision(std::numeric_limits<float>::digits10);    
    for(int k = 0; k < len; k++) {
        std::cout << scientific << ref[re(k)] << "    "<< ref[im(k)] << std::endl;
    }    
#endif
    //statistics
    
    double error  = 0.0;
    double error_ = 0.0;
    int idx = 0;
    for(int i = 0; i < len*2; i++) {
        if(fabs(out[i] - ref[i]) > error) {
            error = fabs(out[i] - ref[i]);
            error_ = ref[i];
            idx = i;
        }
    }
    cout.precision(std::numeric_limits<double>::digits10);
    cout << scientific << "FFTL/FFTW Backward Float:" << error/error_ << "@" << idx << endl;    
    
    
    // clear the site
    
    fftwf_destroy_plan(pFFTW);
    fftwf_free(inFFTW);
    fftwf_free(outFFTW);
    
#endif   
 
    delete[] in;
    delete[] out;
    delete[] ref;    
}




void test_chirpz_d(int len) {
    
    double * in = new double[len*2];
    double * out = new double[len*2];
    double * ref = new double[len*2];    
    
    srand(time(NULL));
    for(int i = 0; i < len*2; i++) {
        in[i] = (double)rand() / (double)RAND_MAX - 0.5;        
    }
    
    unique_ptr<CFastBluestein> pChirpZ(new CFastBluestein);
    pChirpZ->chirpz_d(len, in, out);
    

#ifdef SHOW_VALUES   
    std::cout << "FFTL:" << std::endl;
    cout.precision(std::numeric_limits<double>::digits10);        
    for(int k = 0; k < len; k++) {
        std::cout << scientific << out[re(k)] << "    "<< out[im(k)] << std::endl;
    }
    std::cout << "FFTW:" << std::endl;    
#endif     
    
    
    
    
    // reference output
#ifdef TEST_USE_FFTW        
    fftw_complex *inFFTW, *outFFTW;
    fftw_plan pFFTW;
    
    inFFTW = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * len);
    outFFTW = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * len);
    pFFTW = fftw_plan_dft_1d(len, inFFTW, outFFTW, FFTW_FORWARD, FFTW_MEASURE);    
    memcpy(inFFTW, in, sizeof(double)*len*Cx2Re);         
    fftw_execute(pFFTW);    
    memcpy(ref, outFFTW, sizeof(double)*len*Cx2Re);
    
#ifdef SHOW_VALUES
    cout.precision(std::numeric_limits<double>::digits10);    
    for(int k = 0; k < len; k++) {
        std::cout << scientific << ref[re(k)] << "    "<< ref[im(k)] << std::endl;
    }    
#endif    
    //statistics
    
    double error  = 0.0;
    double error_ = 0.0;
    int idx = 0;
    for(int i = 0; i < len*2; i++) {
        if(fabs(out[i] - ref[i]) > error) {
            error = fabs(out[i] - ref[i]);
            error_ = ref[i];
            idx = i;
        }
    }
    cout.precision(std::numeric_limits<double>::digits10);
    cout << scientific << "FFTL/FFTW Forward Double:" << error/error_ << "@" << idx << endl;    
    
    
    // clear the site
    
    fftw_destroy_plan(pFFTW);
    fftw_free(inFFTW);
    fftw_free(outFFTW);
#endif    
    
    delete[] in;
    delete[] out;
    delete[] ref;    
}





void test_ichirpz_d(int len) {
    
    double * in = new double[len*2];
    double * out = new double[len*2];
    double * ref = new double[len*2];    
    
    srand(time(NULL));
    for(int i = 0; i < len*2; i++) {
        in[i] = (double)rand() / (double)RAND_MAX - 0.5;        
    }
    
    unique_ptr<CFastBluestein> pChirpZ(new CFastBluestein);
    pChirpZ->ichirpz_d(len, in, out);
    

#ifdef SHOW_VALUES   
    std::cout << "FFTL:" << std::endl;
    cout.precision(std::numeric_limits<double>::digits10);        
    for(int k = 0; k < len; k++) {
        std::cout << scientific << out[re(k)] << "    "<< out[im(k)] << std::endl;
    }
    std::cout << "FFTW:" << std::endl;    
#endif     
    
    
    
    
    // reference output
#ifdef TEST_USE_FFTW        
    fftw_complex *inFFTW, *outFFTW;
    fftw_plan pFFTW;
    
    inFFTW = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * len);
    outFFTW = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * len);
    pFFTW = fftw_plan_dft_1d(len, inFFTW, outFFTW, FFTW_BACKWARD, FFTW_MEASURE);    
    memcpy(inFFTW, in, sizeof(double)*len*Cx2Re);         
    fftw_execute(pFFTW);    
    memcpy(ref, outFFTW, sizeof(double)*len*Cx2Re);
    
#ifdef SHOW_VALUES
    cout.precision(std::numeric_limits<double>::digits10);    
    for(int k = 0; k < len; k++) {
        std::cout << scientific << ref[re(k)] << "    "<< ref[im(k)] << std::endl;
    }    
#endif    
    //statistics
    
    double error  = 0.0;
    double error_ = 0.0;
    int idx = 0;
    for(int i = 0; i < len*2; i++) {
        if(fabs(out[i] - ref[i]) > error) {
            error = fabs(out[i] - ref[i]);
            error_ = ref[i];
            idx = i;
        }
    }
    cout.precision(std::numeric_limits<double>::digits10);
    cout << scientific << "FFTL/FFTW Backward Double:" << error/error_ << "@" << idx << endl;    
    
    
    // clear the site
    
    fftw_destroy_plan(pFFTW);
    fftw_free(inFFTW);
    fftw_free(outFFTW);
#endif    
    
    delete[] in;
    delete[] out;
    delete[] ref;    
}



//FFTL/FFTW Forward Float :-2.363980519633000e-02 @130942
//FFTL/FFTW Forward Double:-6.304069676320861e-10 @118993









/*
 * 
 */
int main(int argc, char** argv) {

    
    test_fft_f(65536);  //pow(2, 25)/2 ~ 256MB memory length of time series
    test_ifft_f(65536);
    test_fft_d(65536);
    test_ifft_d(65536);    
    
    test_chirpz_f(65539);
    test_ichirpz_f(65539);    
    test_chirpz_d(65539);
    test_ichirpz_d(65539);    
    
    return 0;
    
    
    // prepare input samples in complex numbers

    float * inSP = new float[N*2];
    float * OutR2 = new float[N*2];
    float * OutR4 = new float[N*2];
    float * OutFFTW = new float[N*2];    

    double * inDP = new double[N*2];
    double * OutR2d = new double[N*2];
    double * OutR4d = new double[N*2];
    double * OutFFTWd = new double[N*2];  
    
    


    // prepare bit-reverse mapping and twiddle factor tables

    unique_ptr<CSort> pSort(new CSort(N));
    unique_ptr<CCoef> pCoef(new CCoef(N));
    unique_ptr<CCoefd> pCoefd(new CCoefd(N));

    unique_ptr<CSort4> pSort4(new CSort4(N));
    unique_ptr<CCoef4> pCoef4(new CCoef4(N));
    unique_ptr<CCoef4d> pCoef4d(new CCoef4d(N));

    

    // radix 2 fft
    pSort->LutComp();
    
    pCoef->CoefComp();
    pCoef->CoefCompI();
    
    pCoefd->CoefComp();
    pCoefd->CoefCompI();


    // radix 4 fft
    pSort4->LutComp();
    
    pCoef4->CoefComp();
    pCoef4->CoefCompI();

    pCoef4d->CoefComp();
    pCoef4d->CoefCompI();    

    
    // create as many transformations as working threads need

    unique_ptr<CTrans>  pFFT  (new CTrans (*pSort, *pCoef));
    unique_ptr<CTrans>  pFFTd (new CTrans (*pSort, *pCoefd));
    unique_ptr<CTrans4> pFFT4 (new CTrans4(*pSort4, *pCoef4));
    unique_ptr<CTrans4> pFFT4d (new CTrans4(*pSort4, *pCoef4d));    
    

  

   
    // prepare time scores
    timespec TimeZero, TimeOne;
    long Elapsed;        
    float * x = NULL;
    double * xDP = NULL;    
    
    

    
    
    fftwf_complex *in, *out;
    fftwf_plan p;
    
    in = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N);
    out = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N);
    p = fftwf_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_MEASURE);
    
    
    
    
    fftw_complex *inDb, *outDb;
    fftw_plan pDb;
    
    inDb = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    outDb = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    pDb = fftw_plan_dft_1d(N, inDb, outDb, FFTW_FORWARD, FFTW_MEASURE);    
    
    
    
    
    for(int j = 0; j < 1; j++) {    
    
    
        srand(time(NULL));
        for(int i = 0; i < N*2; i++) {
            inSP[i] = (float) ((double)rand() / (double)RAND_MAX - 0.5);        
        }
        srand(time(NULL));
        for(int i = 0; i < N*2; i++) {
            inDP[i] = (double)rand() / (double)RAND_MAX - 0.5;        
        }


    
    


    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeZero);    
    for(int i = 0; i < N_TEST_RUN; i++) {
        memcpy(pFFT->GetInputFloat(), inSP, sizeof(float)*N*2);        
        pFFT->InverseFloat();
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeOne);
    Elapsed = (TimeOne.tv_sec - TimeZero.tv_sec) * 1000000000 + (TimeOne.tv_nsec - TimeZero.tv_nsec);
    memcpy(OutR2, pFFT->GetOutputFloat(), sizeof(float)*N*2);
    
    // around -140dB precision...
    x = pFFT->GetOutputFloat();
#ifdef SHOW_VALUES
    cout.precision(std::numeric_limits<float>::digits10);
    for(int i = 0; i < N; i++) {
        cout << scientific << x[i*2+0] << " " << x[i*2+1] << endl;        
    }
#endif
    cout << "R2 SP Time Elapsed:" << Elapsed << " ns" << endl; 

    
    
    
    


    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeZero);
    for(int i = 0; i < N_TEST_RUN; i++) {
        memcpy(pFFT4->GetInputFloat(), inSP, sizeof(float)*N*2);        
        pFFT4->InverseFloat();
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeOne);
    Elapsed = (TimeOne.tv_sec - TimeZero.tv_sec) * 1000000000 + (TimeOne.tv_nsec - TimeZero.tv_nsec);
    memcpy(OutR4, pFFT4->GetOutputFloat(), sizeof(float)*N*2);
    
    
    // around -140dB precision...
    x = pFFT4->GetOutputFloat();
#ifdef SHOW_VALUES
    cout.precision(std::numeric_limits<float>::digits10);
    for(int i = 0; i < N; i++) {
        cout << scientific << x[i*2+0] << " " << x[i*2+1] << endl;        
    }
#endif
    cout << "R4 SP Time Elapsed:" << Elapsed << " ns" << endl; 

    
    


    


    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeZero);
    
    for(int i = 0; i < N_TEST_RUN; i++) {
        memcpy(pFFTd->GetInputDouble(), inDP, sizeof(double)*N*2);        
        pFFTd->InverseDouble();
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeOne);
    Elapsed = (TimeOne.tv_sec - TimeZero.tv_sec) * 1000000000 + (TimeOne.tv_nsec - TimeZero.tv_nsec);
    memcpy(OutR2d, pFFTd->GetOutputDouble(), sizeof(double)*N*2);    
    
    //around -320dB precision...
    xDP = pFFTd->GetOutputDouble();
#ifdef SHOW_VALUES    
    cout.precision(std::numeric_limits<double>::digits10);
    for(int i = 0; i < N; i++) {
        cout << scientific << xDP[i*2+0] << " " << xDP[i*2+1] << endl;        
    }
#endif
    cout << "R2 DP Time Elapsed:" << Elapsed << " ns" << endl; 
    
    
    


    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeZero);    
    for(int i = 0; i < N_TEST_RUN; i++) {
        memcpy(pFFT4d->GetInputDouble(), inDP, sizeof(double)*N*2);        
        pFFT4d->InverseDouble();
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeOne);
    Elapsed = (TimeOne.tv_sec - TimeZero.tv_sec) * 1000000000 + (TimeOne.tv_nsec - TimeZero.tv_nsec);
    memcpy(OutR4d, pFFT4d->GetOutputDouble(), sizeof(double)*N*2);        
    
    //around -320dB precision...
    xDP = pFFT4d->GetOutputDouble();
#ifdef SHOW_VALUES    
    cout.precision(std::numeric_limits<double>::digits10);
    for(int i = 0; i < N; i++) {
        cout << scientific << xDP[i*2+0] << " " << xDP[i*2+1] << endl;        
    }
#endif
    cout << "R4 DP Time Elapsed:" << Elapsed << " ns" << endl; 

    
    
    

    
    
    
   

    

    
    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeZero);
    for(int i = 0; i < N_TEST_RUN; i++) {
        memcpy(in, inSP, sizeof(float)*N*Cx2Re);         
        fftwf_execute(p);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeOne);
    Elapsed = (TimeOne.tv_sec - TimeZero.tv_sec) * 1000000000 + (TimeOne.tv_nsec - TimeZero.tv_nsec);
    memcpy(OutFFTW, out, sizeof(float)*N*2);
    
    cout << "FFTWf Time Elapsed:" << Elapsed << " ns" << endl;     
    
    
 
    


    
    
    
    

    

    
    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeZero);
    for(int i = 0; i < N_TEST_RUN; i++) {
        memcpy(inDb, inDP, sizeof(double)*N*Cx2Re);         
        fftw_execute(pDb);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeOne);
    Elapsed = (TimeOne.tv_sec - TimeZero.tv_sec) * 1000000000 + (TimeOne.tv_nsec - TimeZero.tv_nsec);
    memcpy(OutFFTWd, outDb, sizeof(double)*N*2);
    
    cout << "FFTW Time Elapsed:" << Elapsed << " ns" << endl;     
    
    

    
    
    
    
    
    
 
    //statistics
    float errR2 = 0.0f;
    float errR4 = 0.0f;
    for(int i = 0; i < N*2; i++) {
        if(fabs(OutR2[i] - OutFFTW[i]) > errR2) {
            errR2 = (OutR2[i] - OutFFTW[i]) / OutFFTW[i];
        }
        if(fabs(OutR4[i] - OutFFTW[i]) > errR4) {
            errR4 = (OutR4[i] - OutFFTW[i]) / OutFFTW[i];
        }        
    }
    cout.precision(std::numeric_limits<float>::digits10);
    cout << scientific << "err R2:" << errR2 << "    "<< "err R4:" << errR4 << endl;
    

    double errR2d = 0.0;
    double errR4d = 0.0;
    for(int i = 0; i < N*2; i++) {
        if(fabs(OutR2d[i] - OutFFTWd[i]) > errR2d) {
            errR2d = (OutR2d[i] - OutFFTWd[i]) / OutFFTWd[i];
        }
        if(fabs(OutR4d[i] - OutFFTWd[i]) > errR4d) {
            errR4d = (OutR4d[i] - OutFFTWd[i]) / OutFFTWd[i];
        }        
    }
    cout.precision(std::numeric_limits<double>::digits10);
    cout << scientific << "err R2:" << errR2d << "    "<< "err R4:" << errR4d << endl;

    }
    
    
    fftwf_destroy_plan(p);
    fftwf_free(in);
    fftwf_free(out);
    
    
    
    fftw_destroy_plan(pDb);
    fftw_free(inDb);
    fftw_free(outDb);    
    
    
    delete[] inSP;    
    delete[] inDP;
    
    delete[] OutR2;
    delete[] OutR4;
    delete[] OutFFTW;

    delete[] OutR2d;
    delete[] OutR4d;
    delete[] OutFFTWd;
    
    return 0;
}



//==5845== Memcheck, a memory error detector
//==5845== Copyright (C) 2002-2011, and GNU GPL'd, by Julian Seward et al.
//==5845== Using Valgrind-3.7.0 and LibVEX; rerun with -h for copyright info
//==5845== Command: ./fftl
//==5845== 
//==5845== 
//==5845== HEAP SUMMARY:
//==5845==     in use at exit: 0 bytes in 0 blocks
//==5845==   total heap usage: 184 allocs, 184 frees, 158,611,480 bytes allocated
//==5845== 
//==5845== All heap blocks were freed -- no leaks are possible
//==5845== 
//==5845== For counts of detected and suppressed errors, rerun with: -v
//==5845== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 4 from 4)