/* 
 * File:   CSort4.h
 * Author: lixun
 *
 * Created on 2015年3月14日, 下午11:25
 */

#ifndef CSORT4_H
#define	CSORT4_H

class CSort4 {
    
public:
    CSort4(int Length);
    virtual ~CSort4();
    
    void LutComp(void);
    void LutReal(void);
    const int * GetLutComp(void) const;
    const int * GetLutReal(void) const;
          int   GetPoints(void) const;
          
private:
    CSort4(const CSort4& orig);
    void TreeQuad(void);

protected:
    int * mTreeQuad;
    int * mLutComp;
    int * mLutReal;
    int   mN;
};

#endif	/* CSORT4_H */

