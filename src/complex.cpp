#include "complex.h"
#include <qmath.h>
#include <tgmath.h>

Complex::Complex() {
        re = 0 ;
        im = 0;
}

Complex::Complex(long double r, long double i)    {
        re = r ;
        im = i ;
}


//double Complex::real(void) {return .re;};
//double Complex::imag(void) {return im;};

long double Complex::cabs(void) {
        return sqrtl(re * re + im * im);
};

long double Complex::mag(void) {
        return cabs();
}

long double Complex::arg(void) {
        return atan2l(im, re);
}

// On entre en coord polaire resultat cartesien
Complex   Complex::cPolaire(void) {
        Complex y;
        y.re = re * cosl( im);
        y.im = re * sinl(im);
        return y;
}


