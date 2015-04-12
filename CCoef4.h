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
    
    CCoef4();
    CCoef4(int Length);
    virtual ~CCoef4();
    
    virtual void CoefComp(void);
    virtual void CoefCompI(void);

    const float * GetCoefComp(int Quadrant) const;
    const float * GetCoefCompI(int Quadrant) const;
    
    int GetPoints(void) const;    
    
    
private:

    float * CoWnHeap1;
    float * CoWnHeap2;
    float * CoWnHeap3;    

    float * CoWnHeap1I;
    float * CoWnHeap2I;
    float * CoWnHeap3I; 
    
    CCoef4(const CCoef4& orig);
    void CoefGenHeap(float * CoWnHeap, const long double * WnQ);
    
    
protected:
        
    int mN;
    int mQuadPts;  //number of points per quadrant. mN/4
    int mHeapPts;  //number of points for heaps of all quadrants    

    
    void CoefGenWnQ1(long double * WnQ1);
    void CoefGenWnQ2(long double * WnQ2, const long double * WnQ1);    
    void CoefGenWnQ3(long double * WnQ3, const long double * WnQ1);

    void CoefGenWnQ1I(long double * WnQ1);
    void CoefGenWnQ2I(long double * WnQ2, const long double * WnQ1);    
    void CoefGenWnQ3I(long double * WnQ3, const long double * WnQ1);
    
};

#endif	/* CCOEFRADIX4_H */

