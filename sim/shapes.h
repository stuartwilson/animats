#include "phyllotaxis.h"

class torus: public phyllotaxis {

public:

    torus(int density){
        initPhyllo(0.,2.*M_PI,density);
    };

    // FOR SPHERE (BUT SET T1 = 2.*M_PI)
    vector<double> func(double x){
        vector<double> f(4,0.);
        double a = 0.25;
        double c = 0.5;
        f[0] =  sin(x)*a+c;             //xF
        f[1] =  cos(x)*a;               //zF
        f[2] =  cos(x)*c;               //xdF
        f[3] =  -a*sin(x);              //zdF
        return f;
    };

};


class sphere: public phyllotaxis {

public:
    sphere(int density){
        initPhyllo(0.,M_PI,density);
    };

    vector<double> func(double x){
        vector<double> f(4,0.);
        f[0] =  sin(x);                 //xF
        f[1] =  cos(x);                 //zF
        f[2] =  cos(x);                 //xdF
        f[3] =  -sin(x);                //zdF
        return f;
    }

};
