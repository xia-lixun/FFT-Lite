/* 
 * File:   CFastBluestein.cpp
 * Author: lixun
 * 
 * Created on 2015年4月9日, 上午12:49
 */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <memory>

#include "CChirpZ.h"
#include "CChirpZd.h"
#include "CFastBluestein.h"


CFastBluestein::CFastBluestein() {
}





void CFastBluestein::chirpz_f(int Points, float * In, float * Out) {

    std::unique_ptr<CChirpZ> pChirpZ_f(new CChirpZ(Points)); 
    
    std::memcpy(pChirpZ_f->GetInputFloat(), In, sizeof(float)*Points*2);
    pChirpZ_f->ConvertFloat(M_PIl);
    std::memcpy(Out, pChirpZ_f->GetOutputFloat(), sizeof(float)*Points*2); 
}




void CFastBluestein::ichirpz_f(int Points, float * In, float * Out) {

    std::unique_ptr<CChirpZ> pChirpZ_f(new CChirpZ(Points)); 
    
    std::memcpy(pChirpZ_f->GetInputFloat(), In, sizeof(float)*Points*2);
    pChirpZ_f->ConvertFloat(-M_PIl);
    std::memcpy(Out, pChirpZ_f->GetOutputFloat(), sizeof(float)*Points*2);     
}




void CFastBluestein::chirpz_d(int Points, double * In, double * Out) {
    
    std::unique_ptr<CChirpZd> pChirpZ_d(new CChirpZd(Points));
    
    std::memcpy(pChirpZ_d->GetInputDouble(), In, sizeof(double)*Points*2);
    pChirpZ_d->ConvertDouble(M_PIl);
    std::memcpy(Out, pChirpZ_d->GetOutputDouble(), sizeof(double)*Points*2);      
}




void CFastBluestein::ichirpz_d(int Points, double * In, double * Out) {

    std::unique_ptr<CChirpZd> pChirpZ_d(new CChirpZd(Points));
    
    std::memcpy(pChirpZ_d->GetInputDouble(), In, sizeof(double)*Points*2);
    pChirpZ_d->ConvertDouble(-M_PIl);
    std::memcpy(Out, pChirpZ_d->GetOutputDouble(), sizeof(double)*Points*2);          
}





CFastBluestein::CFastBluestein(const CFastBluestein& orig) {
}

CFastBluestein::~CFastBluestein() {
}

