/* 
 * File:   CCoefRadix4.h
 * Author: lixun
 *
 * Created on 2015年3月10日, 下午7:53
 */

#ifndef CCOEFRADIX4_H
#define	CCOEFRADIX4_H



class CCoef4 {

public:
    CCoef4(int Length);
    virtual ~CCoef4();
    
    void CoefComp(void);
    void CoefCompI(void);

    const float * GetCoefComp(int Quadrant) const;
    const float * GetCoefCompI(int Quadrant) const;
    int           GetPoints(void) const;    
    
private:
    CCoef4(const CCoef4& orig);
    void CoefGen(int Quadrant, float * CoWnHeap, double Math2);    
    
protected:
    double * Wn;
    
    float * CoWnHeap1;
    float * CoWnHeap2;
    float * CoWnHeap3;    

    float * CoWnHeap1I;
    float * CoWnHeap2I;
    float * CoWnHeap3I;    
    
    int     mN;
    int     mQuadPts;  //number of points per quadrant. mN/4
    int     mHeapPts;  //number of points for heaps of all quadrants    
 
};

#endif	/* CCOEFRADIX4_H */

