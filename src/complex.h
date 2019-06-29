#ifndef COMPLEX_H
#define COMPLEX_H

#include <qmath.h>

class Complex {
        public:
                long double re, im;

        public:
                Complex(long double r, long double i) ;
                Complex(void);


                Complex&   operator += (const Complex& rhs) {
                        re += rhs.re ;
                        im += rhs.im;
                        return  *this;
                }
                Complex  operator +( const Complex& rhs) {
                        Complex cRet;
                        cRet.re = rhs.re + this->re ;
                        cRet.im = rhs.im + this->im;
                        return  cRet;
                }
                Complex& operator -= (const Complex& rhs) {
                        re -= rhs.re;
                        im -= rhs.im;
                        return *this;
                }
                Complex operator - (const Complex& rhs) {
                          Complex cRet;
                          cRet.re = -rhs.re + this->re ;
                          cRet.im = -rhs.im + this->im;
                          return  cRet;
                }
                Complex& operator *= (const  Complex& rhs) {
                        long double tmp = re * rhs.re - im * rhs.im;
                        im = im * rhs.re + re * rhs.im;
                        re = tmp;
                        return *this;
                }
         /*      Complex& operator *= (const  long double& rhs ) { marche pas
                   Complex cRet;
                        cRet.im = im * rhs;
                      cRet.re=  re * rhs;
                        return cRet;
                }*/
                Complex operator * (const  Complex& rhs) {
                         Complex cRet;
                         cRet.im = im * rhs.re + re * rhs.im;
                        cRet.re =  re * rhs.re - im * rhs.im;
                        return cRet;
                }
                Complex& operator /= (const  Complex& rhs) {
                        long double dr, di;
                        long double den = rhs.re * rhs.re + rhs.im * rhs.im;
                        if (den == 0) {
                                re = 0.0;
                                im = 0.0;
                        } else {
                                dr = re;
                                di = im;
                                re = (dr * rhs.re + di * rhs.im) / den;
                                im = ( di * rhs.re - dr * rhs.im) / den;
                        }
                        return *this;
                }
                Complex operator / (const  Complex& rhs) {
                        Complex cRet;
                        long double den = rhs.re * rhs.re + rhs.im * rhs.im;
                        if (den == 0) {
                                cRet = Complex(0, 0);
                        } else {
                                cRet.re = (re * rhs.re + im * rhs.im) / den;
                                cRet.im = ( im * rhs.re - re * rhs.im) / den;
                        }
                        return cRet;
                }


                //  Complex & operator= (const  double& r,const  double& i) {re=r; im=i; return *this;}
                long  double cabs(void);
                long double mag(void);
                long double arg(void);
                Complex  cPolaire(void);
};

#endif // COMPLEX_H
