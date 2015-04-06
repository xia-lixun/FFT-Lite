/* 
 * File:   CFastFourier.cpp
 * Author: lixun
 * 
 * Created on 2015年4月1日, 上午12:09
 */

#include <cstdlib>
#include <cstring>
#include <memory>

#include "CSort.h"
#include "CSort4.h"
#include "CCoef.h"
#include "CCoef4.h"
#include "CTrans.h"
#include "CTrans4.h"

#include "CFastFourier.h"





bool CFastFourier::IsRadix4(int Length) {
    
    bool isRadix4 = false;
    while(Length > 0) {
        Length >>= 2;
        if(Length == 1) {
            isRadix4 = true;
            break;
        }
    }
    return isRadix4;    
}






void CFastFourier::fft_f(int Points, float * In, float * Out) {
    
    if(IsRadix4(Points)) {
        
        std::unique_ptr<CSort4> pSort4(new CSort4(Points));
        std::unique_ptr<CCoef4> pCoef4(new CCoef4(Points)); 
        
        pSort4->LutComp();
        pCoef4->CoefComp();

        std::unique_ptr<CTrans4> pFFT4 (new CTrans4(*pSort4, *pCoef4));        
        
        
        std::memcpy(pFFT4->GetInputFloat(), In, sizeof(float)*Points*2);        
        pFFT4->ForwardFloat();       
        std::memcpy(Out, pFFT4->GetOutputFloat(), sizeof(float)*Points*2);
    }
    else {
        
        std::unique_ptr<CSort> pSort(new CSort(Points));
        std::unique_ptr<CCoef> pCoef(new CCoef(Points));
        
        pSort->LutComp();   
        pCoef->CoefComp();

        std::unique_ptr<CTrans>  pFFT  (new CTrans (*pSort, *pCoef));        
        
        
        std::memcpy(pFFT->GetInputFloat(), In, sizeof(float)*Points*2);        
        pFFT->ForwardFloat();
        std::memcpy(Out, pFFT->GetOutputFloat(), sizeof(float)*Points*2);
    }
}





void CFastFourier::ifft_f(int Points, float * In, float * Out) {
    
    if(IsRadix4(Points)) {
        
        std::unique_ptr<CSort4> pSort4(new CSort4(Points));
        std::unique_ptr<CCoef4> pCoef4(new CCoef4(Points)); 
        
        pSort4->LutComp();
        pCoef4->CoefCompI();

        std::unique_ptr<CTrans4> pIFFT4 (new CTrans4(*pSort4, *pCoef4));        
        
        
        std::memcpy(pIFFT4->GetInputFloat(), In, sizeof(float)*Points*2);        
        pIFFT4->InverseFloat();       
        std::memcpy(Out, pIFFT4->GetOutputFloat(), sizeof(float)*Points*2);
    }
    else {
        
        std::unique_ptr<CSort> pSort(new CSort(Points));
        std::unique_ptr<CCoef> pCoef(new CCoef(Points));
        
        pSort->LutComp();
        pCoef->CoefCompI();

        std::unique_ptr<CTrans>  pIFFT  (new CTrans (*pSort, *pCoef));        
        
        
        std::memcpy(pIFFT->GetInputFloat(), In, sizeof(float)*Points*2);        
        pIFFT->InverseFloat();
        std::memcpy(Out, pIFFT->GetOutputFloat(), sizeof(float)*Points*2);
    }
}





void CFastFourier::fft_d(int Points, double * In, double * Out) {
    
    if(IsRadix4(Points)) {
        
        std::unique_ptr<CSort4> pSort4(new CSort4(Points));
        std::unique_ptr<CCoef4d> pCoef4d(new CCoef4d(Points));       
        
        pSort4->LutComp();
        pCoef4d->CoefComp();

        std::unique_ptr<CTrans4> pFFT4d (new CTrans4(*pSort4, *pCoef4d));    
        
        
        std::memcpy(pFFT4d->GetInputDouble(), In, sizeof(double)*Points*2);        
        pFFT4d->ForwardDouble();       
        std::memcpy(Out, pFFT4d->GetOutputDouble(), sizeof(double)*Points*2);
    }
    else {
        
        std::unique_ptr<CSort> pSort(new CSort(Points));
        std::unique_ptr<CCoefd> pCoefd(new CCoefd(Points));       
        
        pSort->LutComp();   
        pCoefd->CoefComp();

        std::unique_ptr<CTrans>  pFFTd  (new CTrans (*pSort, *pCoefd));        
        
        
        std::memcpy(pFFTd->GetInputDouble(), In, sizeof(double)*Points*2);        
        pFFTd->ForwardDouble();
        std::memcpy(Out, pFFTd->GetOutputDouble(), sizeof(double)*Points*2);
    }
}


void CFastFourier::ifft_d(int Points, double * In, double * Out) {
    
    if(IsRadix4(Points)) {
        
        std::unique_ptr<CSort4> pSort4(new CSort4(Points));
        std::unique_ptr<CCoef4d> pCoef4d(new CCoef4d(Points));       
        
        pSort4->LutComp();
        pCoef4d->CoefCompI();

        std::unique_ptr<CTrans4> pIFFT4d (new CTrans4(*pSort4, *pCoef4d));    
        
        
        std::memcpy(pIFFT4d->GetInputDouble(), In, sizeof(double)*Points*2);        
        pIFFT4d->InverseDouble();       
        std::memcpy(Out, pIFFT4d->GetOutputDouble(), sizeof(double)*Points*2);
    }
    else {
        
        std::unique_ptr<CSort> pSort(new CSort(Points));
        std::unique_ptr<CCoefd> pCoefd(new CCoefd(Points));       
        
        pSort->LutComp();   
        pCoefd->CoefCompI();

        std::unique_ptr<CTrans>  pIFFTd  (new CTrans (*pSort, *pCoefd));        
        
        
        std::memcpy(pIFFTd->GetInputDouble(), In, sizeof(double)*Points*2);        
        pIFFTd->InverseDouble();
        std::memcpy(Out, pIFFTd->GetOutputDouble(), sizeof(double)*Points*2);
    }
}






CFastFourier::CFastFourier() {
}

CFastFourier::CFastFourier(const CFastFourier& orig) {
}

CFastFourier::~CFastFourier() {
}


