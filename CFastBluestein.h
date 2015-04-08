/* 
 * File:   CFastBluestein.h
 * Author: lixun
 *
 * Created on 2015年4月9日, 上午12:49
 */

#ifndef CFASTBLUESTEIN_H
#define	CFASTBLUESTEIN_H





class CFastBluestein {
    
public:
    CFastBluestein();
    virtual ~CFastBluestein();
    
    void chirpz_f(int Points, float * In, float * Out);
    void ichirpz_f(int Points, float * In, float * Out);
    void chirpz_d(int Points, double * In, double * Out);
    void ichirpz_d(int Points, double * In, double * Out);    

    
private:
    CFastBluestein(const CFastBluestein& orig);


protected:    

};

#endif	/* CFASTBLUESTEIN_H */

