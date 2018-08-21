#include <math.h>
#include <vector>

using namespace std;

class sphere {

public:

    vector<vector<int> > sphereC, sphereM;
    vector<vector<double> > sphereX, sphereS;
    int sphereN;
    double radius;

    sphere(void){

    }

    sphere(int n){

        initSphere(n,1.0);
    }

    void initSphere(int n, double radius){

        /*
         cube=sphere indexing system based on
         http://www.iquilezles.org/www/articles/patchedsphere/patchedsphere.htm
         and spherical coordinates based on
         http://mathworld.wolfram.com/SphericalCoordinates.html
         */
        this->sphereN = n;
        this->radius = radius;

        int N = n*n*6;
        vector<double> px(3,0.);
        vector<double> ps(2,0.);
        vector<int> pc(3,0);
        for(int f=0;f<6;f++){
            for(int s=0;s<n;s++){
                double sd=((double)s+0.5)/(double)n;
                for(int t=0;t<n;t++){
                    double td=((double)t+0.5)/(double)n;
                    double x=1.-2.*sd;
                    double y=1.-2.*td;
                    double k=1./sqrt(x*x+y*y+1);
                    x*=k;
                    y*=k;
                    switch(f){
                        case(0):{ //+z
                            px[0]=-x;
                            px[1]=-y;
                            px[2]=k;
                        } break;
                        case(1):{ //-y
                            px[0]=-x;
                            px[1]=-k;
                            px[2]=-y;
                        } break;
                        case(2):{ //-z
                            px[0]=-x;
                            px[1]=y;
                            px[2]=-k;
                        } break;
                        case(3):{ //+y
                            px[0]=-x;
                            px[1]=k;
                            px[2]=y;
                        } break;
                        case(4):{ //-x
                            px[0]=-k;
                            px[1]=x;
                            px[2]=-y;
                        } break;
                        case(5):{ //+x
                            px[0]=k;
                            px[1]=-x;
                            px[2]=-y;
                        } break;
                    }

                    px[0]*=radius;
                    px[1]*=radius;
                    px[2]*=radius;

                    // Cartesian
                    sphereX.push_back(px);

                    // Spherical

                    // theta (round equator) in 0->2pi
                    ps[0]=atan2(px[1],px[0])+M_PI;
                    // CONSIDER WRAPPING ANGLES HERE a-2pi*floor(a/(2pi))
                    // phi (azimuthal) from +z in 0->pi
                    ps[1]=acos(px[2]);
                    sphereS.push_back(ps);

                    // Cube
                    pc[0]=f;
                    pc[1]=s;
                    pc[2]=t;
                    sphereC.push_back(pc);
                }
            }
        }
        sphereM.resize(N);
        for(int i=0;i<N;i++){
            sphereM[i].resize(4);
        }

        int n2 = n*n;

        for(int i=0;i<6;i++){
            for(int s=0;s<n-1;s++){
                for(int t=0;t<n;t++){
                    sphereM[i*n2+s*n+t][0]=i*n2+(s+1)*n+t; //+x
                }
            }
            for(int s=1;s<n;s++){
                for(int t=0;t<n;t++){
                    sphereM[i*n2+s*n+t][2]=i*n2+(s-1)*n+t; //-x
                }
            }
        }

        for(int i=0;i<6;i++){
            for(int s=0;s<n;s++){
                for(int t=0;t<n-1;t++){
                    sphereM[i*n2+s*n+t][1]=i*n2+s*n+t+1; //+y
                }
                for(int t=1;t<n;t++){
                    sphereM[i*n2+s*n+t][3]=i*n2+s*n+t-1; //-y
                }

            }
        }

        int k = n-1;
        for(int i=0;i<n;i++){
            int j = n-i-1;
            // +x
            sphereM[0*n2+n*k+i][0]=  5*n2+n*i+k;
            sphereM[1*n2+n*k+i][0]=  5*n2+n*0+i;
            sphereM[2*n2+n*k+i][0]=  5*n2+n*j+0;
            sphereM[3*n2+n*k+i][0]=  5*n2+n*k+j;
            sphereM[4*n2+n*k+i][0]=  1*n2+n*0+i;
            sphereM[5*n2+n*k+i][0]=  3*n2+n*k+j;

            // +y
            sphereM[0*n2+n*i+k][1]=  3*n2+n*i+0;
            sphereM[1*n2+n*i+k][1]=  0*n2+n*i+0;
            sphereM[2*n2+n*i+k][1]=  1*n2+n*i+0;
            sphereM[3*n2+n*i+k][1]=  2*n2+n*i+0;
            sphereM[4*n2+n*i+k][1]=  0*n2+n*0+j;
            sphereM[5*n2+n*i+k][1]=  0*n2+n*k+i;

            // -x
            sphereM[0*n2+n*0+i][2]=  4*n2+n*j+k;
            sphereM[1*n2+n*0+i][2]=  4*n2+n*k+i;
            sphereM[2*n2+n*0+i][2]=  4*n2+n*i+0;
            sphereM[3*n2+n*0+i][2]=  4*n2+n*0+j;
            sphereM[4*n2+n*0+i][2]=  3*n2+n*0+j;
            sphereM[5*n2+n*0+i][2]=  1*n2+n*k+i;

            // -y
            sphereM[0*n2+n*i+0][3]=  1*n2+n*i+k;
            sphereM[1*n2+n*i+0][3]=  2*n2+n*i+k;
            sphereM[2*n2+n*i+0][3]=  3*n2+n*i+k;
            sphereM[3*n2+n*i+0][3]=  0*n2+n*i+k;
            sphereM[4*n2+n*i+0][3]=  2*n2+n*0+j;
            sphereM[5*n2+n*i+0][3]=  2*n2+n*k+j;
        }

    }




    ~sphere(void){

    };


    /*
    vector<vector<int> > sphereC, sphereM;
    vector<vector<double> > sphereX, sphereS;
    int sphereN;
    double radius;
    sphere(void);
    sphere(int);
    void initSphere(int,double);
    virtual ~sphere();
     */

};


/*
#ifndef ____sphere__
#define ____sphere__

#include <math.h>
#include <vector>

using namespace std;

class sphere {
    
public:
    vector<vector<int> > sphereC, sphereM;
    vector<vector<double> > sphereX, sphereS;
    int sphereN;
    double radius;
    sphere(void);
    sphere(int);
    void initSphere(int,double);
    virtual ~sphere();

};
*/
//#endif /* defined(____sphere__) */
