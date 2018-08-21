///
///  @file shapeMatch.h
///  @brief Shape match algorithm
///  @author Hiroki Urashima, Stuart P Wilson
///  @date 01/04/2013
///

#ifndef ____shapeMatch__
#define ____shapeMatch__

#include <vector>
#include <armadillo>

using std::vector;
using arma::vec;
using arma::mat;

typedef enum{
    Rigid,    // Rigid body
    Linear,   // Linear deformation
    Quadratic // Quadratic deformation
} ShapeMatchType;

class ShapeMatchParameters{
public:
    ShapeMatchParameters();
    virtual ~ShapeMatchParameters();
    double dt;
    vec external;
    double alpha;               // rigidity (from 0  to 1, where 1 is rigid)
    double beta;                // non-linearity (from 0 to 1, where 1 is completely non-linear)
    bool volumeConservation;
    bool allowFlip;
    ShapeMatchType type;
};

class ShapeMatch{
public:
    ShapeMatch();
    virtual ~ShapeMatch();
    void initShape(const vector<vec>,
              const vector<double>,
              const vector<vector<int> >,
              const vector<int>);
    
    void externalForces();
    void externalForces(double,double,double);
    vector<vector<double> > getNormals(void);
    void inflate();
    vector<double> getCoG(void);
    void calculateCollisionDetectionWith(ShapeMatch s);
    void projectPositions();
    void integrate();
    
    void relocate(const vector<vec> );
    void updateBonePosition(const vector<vec> );
    void iterate();
    ShapeMatchParameters params;
    
    int Npoint;              // number of points in cloud
    vector<vec> Xnow;      // current points loci
    vector<vec> Xnew;   // projected point loci
    vector<vec> Xtemp;  // temporary copy of projected loci
    vector<vec> Vel;    // velocities of points
    vector<int> boneIDs;// indices of bones (for external reference only)
    void reset();
    vector<bool> free;
    int Nclusters;
    vector<vector<int> > clusterIDs;    // indices of points per cluster
    vector<vec> Xorg;   // original loci
    vector<vector<vec> > Xgoal;         // goal loci
        void precompute();
private:

    void precompute(int);

    vector<vector<double> > Mass;       // mass points
    
    vector<mat> m_q;    // precomputed
    vector<mat> m_qQ;   // precomputed
    vector<mat> m_Aqq;  // precomputed
    vector<mat> m_AqqQ; // precomputed
    double third;
};

#endif // ____shapeMatch__
