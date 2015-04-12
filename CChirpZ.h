/* 
 * File:   CChirpZ.h
 * Author: lixun
 *
 * Created on 2015年4月5日, 下午2:52
 */

#ifndef CCHIRPZ_H
#define	CCHIRPZ_H

class CChirpZ {
    
public:
    
    CChirpZ(int Points);
    virtual ~CChirpZ();
    
    float * GetInputFloat(void);
    float * GetOutputFloat(void);    
    void ConvertFloat(long double MathPiWithSign);  // M_PIl = forward, -M_PIl = backward        
    
private:
    
    CChirpZ(const CChirpZ& orig);
    
protected:
    
    int mPoints;
    int mPointsChirpZ;  // radix 2 and greater than or equal to 2 * mPoints - 1
    
    float * C;
    float * X;
    
    float * Xs;
    float * Cs;
   
};

#endif	/* CCHIRPZ_H */

