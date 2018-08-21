/*!
 * A template RD system
 */
#include "shapeMatch.h"
#include "sphere.h"

#include "morph/display.h"
#include <iostream>
#include <vector>
#include <string>

using namespace std;



// PLANE BOUNDARY, E.G., FOR FLOOR
void boundary(ShapeMatch& a, double loc, int axis, double bounce, int contactFlag){

    if (loc<0.){
        for(int j=0;j<a.Npoint;j++){
            if (a.Xtemp[j][axis]<loc+a.rad[j]){
                vec normal = zeros<vec>(3);
                normal(axis) = 1.0;
                vec vNew = -2.*dot(a.Vel[j],normal)*normal+a.Vel[j];
                a.Xnew[j]=a.Xnow[j]+vNew*a.params.dt*bounce;
                if(a.Xnew[j][axis]<loc+a.rad[j]){
                    a.Xnew[j][axis]=loc+a.rad[j];
                    //a.contacts[j]=contactFlag;
                }
            }
        }
    }else{
        for(int j=0;j<a.Npoint;j++){
            if (a.Xtemp[j][axis]>loc-a.rad[j]){
                vec normal = zeros<vec>(3);
                normal(axis) = -1.0;
                vec vNew = -2.*dot(a.Vel[j],normal)*normal+a.Vel[j];
                a.Xnew[j]=a.Xnow[j]+vNew*a.params.dt*bounce;
                if(a.Xnew[j][axis]>loc-a.rad[j]){
                    a.Xnew[j][axis]=loc-a.rad[j];
                    //a.contacts[j]=contactFlag;
                }
            }
        }
    }
}




int main (int argc, char **argv)
{
    if (argc < 2) {
        cerr << "\nUsage: ./build/sim/process w0\n\n";
        cerr << "Be sure to run from the base NeoArealize source directory.\n";
        return -1;
    }

    // Set RNG seed
    int rseed = 1;
    srand(rseed);

    // Create some displays
    vector<morph::Gdisplay> displays;
    vector<double> fix(3, 0.0);
    vector<double> eye(3, 0.0);
    vector<double> rot(3, 0.0);

    double rhoInit = 5.;
    string worldName(argv[1]);
    string winTitle = worldName + ": animat arena";
    displays.push_back (morph::Gdisplay (500, 500, 100, 0, winTitle.c_str(), rhoInit, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();
    // Instantiate the model object
    ShapeMatch SM;

    sphere S(7);
    int N = S.sphereN*S.sphereN*6;
    vector<double> mass (N, 1.0/(double)N);
    vector<vector<int> > clusters(1);
    clusters.resize(1);
    clusters[0].resize(N,0);
    for(int i=0;i<N;i++){clusters[0][i]=i;}
    vector<int> boneIDs;
    vector<arma::vec> X;
    for(int i=0;i<N;i++){X.push_back(S.sphereX[i]);};

    try {
        SM.params.alpha = 0.2;
        SM.params.beta = 1.0;
        SM.params.type = Quadratic;
        SM.params.dt = 0.01;
        SM.params.allowFlip = true;
        SM.params.volumeConservation = true;
        SM.initShape (X, mass, clusters, boneIDs, 0.05 );
        SM.precompute();

    } catch (const exception& e) {
        cerr << "Exception initialising SM object: " << e.what() << endl;
    }



    // Start the loop
    bool finished = false;
    while (!finished) {
        // Step the model
        try {
            SM.externalForces(0.,0.,-10.);
            SM.projectPositions();
            boundary(SM, -1.5, 2, 0.5, 0);
            SM.integrate();
            //SM.step();
        } catch (const exception& e) {
            cerr << "Caught exception calling SM.step(): " << e.what() << endl;
            finished = true;
        }

        displays[0].resetDisplay (fix, eye, rot);
        try {
            displays[0].resetDisplay (fix, eye, rot);
            vector<double> cl(3,0.5);

            for(int i=0;i<N;i++){
            displays[0].drawSphere(SM.Xnow[i](0),SM.Xnow[i](1),SM.Xnow[i](2),SM.rad[i],cl,9);
            }

            displays[0].redrawDisplay();
        } catch (const exception& e) {
            cerr << "Caught exception calling SM.plot(): " << e.what() << endl;
            finished = true;
        }
        /*
        try {
            SM.plot (displays);
            // Save some frames
            if (SM.stepCount % 100 == 0) {
                //RD.saveStuff();
            }

        } catch (const exception& e) {
            cerr << "Caught exception calling SM.plot(): " << e.what() << endl;
            finished = true;
        }


        // Halt after how every many iterations suits your model:
        if (SM.stepCount > 100) {
            finished = true;
        }
         */
    }

    return 0;
};
