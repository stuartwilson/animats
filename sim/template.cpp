/*!
 * A template RD system
 */
#include "rd_template.h"

#include "morph/display.h"
#include <iostream>
#include <vector>
#include <string>

using namespace std;

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
    eye[2] = -0.4;
    vector<double> rot(3, 0.0);

    double rhoInit = 1.5;
    string worldName(argv[1]);
    string winTitle = worldName + ": window name";
    displays.push_back (morph::Gdisplay (1020, 300, 100, 0, winTitle.c_str(), rhoInit, 0.0, 0.0));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    // SW - Contours
    winTitle = worldName + ": contours";
    displays.push_back (morph::Gdisplay (500, 500, 100, 900, winTitle.c_str(), rhoInit, 0.0, 0.0, displays[0].win));
    displays.back().resetDisplay (fix, eye, rot);
    displays.back().redrawDisplay();

    // Instantiate the model object
    RD_Template RD;

    // Do any modifications to RD parameters here.
    // RD.myparam = 3.4;

    // Call the init function, which can allocate variables and run
    // any pre-stepping computations.
    try {
        RD.init (displays);
    } catch (const exception& e) {
        cerr << "Exception initialising RD object: " << e.what() << endl;
    }

    // Start the loop
    bool finished = false;
    while (!finished) {
        // Step the model
        try {
            RD.step();
        } catch (const exception& e) {
            cerr << "Caught exception calling RD.step(): " << e.what() << endl;
            finished = true;
        }

        displays[0].resetDisplay (fix, eye, rot);
        try {
            RD.plot (displays);
            // Save some frames
            if (RD.stepCount % 100 == 0) {
                //RD.saveStuff();
            }

        } catch (const exception& e) {
            cerr << "Caught exception calling RD.plot(): " << e.what() << endl;
            finished = true;
        }

        // Halt after how every many iterations suits your model:
        if (RD.stepCount > 100) {
            finished = true;
        }
    }

    return 0;
};
