/* 
 * File:   Config.h
 * Author: lixun
 *
 * Created on January 24, 2015, 5:40 PM
 */

#ifndef CONFIG_H
#define	CONFIG_H




//////////////////////
// user definitions //
//////////////////////
#define  N          (256)              //N must be greater than or equal to 4
//
// make sure that number of samples must be power of two, if
// not, take the minimal power of two that is greater than
// the number of samples as operation length.
//
//if ((NumSamples & (NumSamples-1)) != 0) {
//    N = 1;
//    while(N < NumSamples) {
//        N <<= 1;
//    }
//} else {
//    N = NumSamples;
//}









////////////////////////
// system definitions //
////////////////////////
#define  Cx2Re       (2)                  //Conversion of Index from real to complex domain
#define  re(x)       ((x) * Cx2Re + 0)    //real index of the complex number indexed by x
#define  im(x)       ((x) * Cx2Re + 1)    //real index of the complex number indexed by x
#define  CacheAlign  (64)                 //memory address being multiple of for cache alignment





#endif	/* CONFIG_H */

