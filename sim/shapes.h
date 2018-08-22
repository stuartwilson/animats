using arma::vec;
using arma::mat;

class phyllotaxis {

public:

    int N;
    vector<vector<double> > X;

    /*
     dummy function to be overwritten by specific shapes
     */
    virtual vector<double> func(double){
        return vector<double> (4,0.);
    }

    double integrand(double t){
        vector<double> f = func(t);
        return 2.0*M_PI * f[0] * pow(pow(f[2],2) + pow(f[3],2),0.5);
    }

    double quad(int res,double t0, double t1){
        vec xint2(res);
        vec yint2(res);
        for(int i=0;i<res;i++){
            double x = t0+(t1-t0)*(double)i/(double)res;
            xint2(i)=x;
            yint2(i)=integrand(x);
        }
        mat z = trapz(xint2,yint2);
        return (double) z(0);
    }

    void init(double t0, double t1, int density){

        double totalArea = fabs(quad(500,t0,t1));

        this->N = floor(density*totalArea);

        this->X.resize(N);

        vector<double> interp (N,0.);
        vector<double> T0(N,0.);
        vector<double> T1(N,0.);

        double dt = t1/(double)N;

        for(int i=0; i<N; i++){
            T0[i] = i*dt;
            T1[i] = (i+1)*dt;
        }

        int res = 50;
        double area1 = quad(res,t0,dt);
        for(int i=0; i<N; i++){
            double area0 = area1;
            area1 += quad(res,T0[i],T1[i]);
            double t = totalArea*((double)i/(double)(N-1));     // interpolation sample points
            double y = (area0*(T1[i]-t)+area1*(t-T0[i]))/dt;    // linear interpolation
            double theta = (double)i*2.39996322972865332;       // fibonnacci golden angle
            X[i].resize(3);
            vector<double> f = func(y);
            X[i][0] = f[0] * cos(theta);
            X[i][1] = f[0] * sin(theta);
            X[i][2] = f[1];
        }
    };

    ~phyllotaxis(void){
    };
};



class torus: public phyllotaxis {

public:

    double a = 0.25;
    double c = 0.5;

    torus(int density){
        init(0.,2.*M_PI,density);
    };

    torus(int density, double a, double c){
        this->a = a;
        this->c = c;
        init(0.,2.*M_PI,density);
    };

    // FOR SPHERE (BUT SET T1 = 2.*M_PI)
    vector<double> func(double x){
        vector<double> f(4,0.);
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
        init(0.,M_PI,density);
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
