#include <vector>
#include <armadillo>

using namespace arma;
using namespace std;

using arma::vec;
using arma::mat;

class phyllotaxis {

public:

    vector<vector<double> > X;
    int N;

    double gamma = 2.39996322972865332; // Golden angle

    phyllotaxis(void){
        initPhyllo(0.,2.*M_PI,50);
    };

    // FOR SPHERE (BUT SET T1 = 2.*M_PI)
    vector<double> func(double x){
        vector<double> f(4,0.);
        double a = 0.25;
        double c = 0.5;
        f[0] =  sin(x)*a+c;          //xF
        f[1] =  cos(x)*a;            //zF
        f[2] =  cos(x)*c;            //xdF
        f[3] =  -a*sin(x);           //zdF
        return f;
    }

    /*
    // FOR SPHERE (BUT SET T1 = M_PI)
    vector<double> func(double x){
        vector<double> f(4,0.);
        f[0] =  sin(x);//*a+c;          //xF
        f[1] =  cos(x);//*a;            //zF
        f[2] =  cos(x);//*c;            //xdF
        f[3] =  -sin(x);           //zdF
        return f;
    }
     */

    double integrand(double t){
        vector<double> f = func(t);
        return 2.0*M_PI * f[0] * pow(pow(f[2],2) + pow(f[3],2),0.5);
    }

    double quad(int res,double t0, double t1){

        vec xint2(res);
        vec yint2(res);
        for(int i=0;i<res;i++){
            double x = t0+(t1-t0)*(double)i/(double)res;
            double y = integrand(x);
            xint2(i)=x;
            yint2(i)=y;
        }
        mat z = trapz(xint2,yint2);
        return (double) z(0);
    }

    void initPhyllo(double t0, double t1, int density){


        double totalArea = fabs(quad(500,t0,t1));

        //Number of points = density * total area rounded down
        this->N = floor(density*totalArea);

        // Equivalent to area function
        vector<double> area(N,0.);
        double inc = t1/(double)N;
        double t1n = inc;

        int res = 50;
        area[0] = quad(res,t0,t1n);
        for(int i=1; i<N; i++){
            double t0n = t1n;
            t1n = t1n+inc;
            //area[i] = quad(res,t0n,t1n) + area[i-1];
            area[i] = quad(res,t0n,t1n) + area[i-1];
        }


    // Equivalent to areaToTFunction in wolfram alpha blog
    // Require to know points at totalarea * i/n


        /*
        t = np.linspace(t0,t1,n);
    areaInterp = np.zeros([n,1]);
    for i in range(n):
    areaInterp[i] =  totalArea*(i/n);
         */

        //X.resize(N);
        vector<double> x(3,0.);
        for(int i=0;i<N;i++){
            double theta = (double)i*gamma; //Iteration * golden angle
            vector<double> f = func(area[i]);
            x[0] = f[0] * cos(theta);
            x[1] = f[0] * sin(theta);
            x[2] = f[1];
            X.push_back(x);
        }
    };

    ~phyllotaxis(void){

    };

};


