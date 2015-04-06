/* 
 * File:   CChirpZd.h
 * Author: lixun
 *
 * Created on 2015年4月5日, 下午6:40
 */

#ifndef CCHIRPZD_H
#define	CCHIRPZD_H

class CChirpZd {
    
public:
    CChirpZd(int Points);
    virtual ~CChirpZd();
    
    double * GetInputDouble(void);
    double * GetOutputDouble(void);    
    void ConvertDouble(void);        
    
    
private:
    CChirpZd(const CChirpZd& orig);
    
protected:
    int mPoints;
    int mPointsChirpZ;  // radix 2 and greater than or equal to 2 * mPoints - 1
    
    double * C;
    double * X;
    
    double * Xs;
    double * Cs;    
};

#endif	/* CCHIRPZD_H */

