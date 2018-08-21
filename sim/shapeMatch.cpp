#include "shapeMatch.h"

using namespace arma;

ShapeMatchParameters::ShapeMatchParameters():

// Default parameters
dt(0.01),
alpha(0.20),
beta(0.20),
type(Quadratic),
volumeConservation(true),
allowFlip(true),
external(){
    external<<0.0<<0.0<<-10.00;
}

ShapeMatchParameters::~ShapeMatchParameters()
{
}

ShapeMatch::ShapeMatch():
Npoint(0),
Xorg(),
Xnow(),
Xnew(),
Xtemp(),
Xgoal(),
clusterIDs(),
boneIDs(),
Mass(),
Vel(),
m_q(),
m_qQ(),
m_Aqq(),
m_AqqQ()
{
}

ShapeMatch::~ShapeMatch()
{
    Npoint=0;
    Xorg.clear();
    Xnow.clear();
    Xnew.clear();
    Xtemp.clear();
    Xgoal.clear();
    clusterIDs.clear();
    boneIDs.clear();
    Mass.clear();
    Vel.clear();
    m_q.clear();
    m_qQ.clear();
    m_Aqq.clear();
    m_AqqQ.clear();
    
}

/// Resets all the member variables
void ShapeMatch::reset(){
    Npoint=0;
    Xorg.clear();
    Xnow.clear();
    Xnew.clear();
    Xtemp.clear();
    Xgoal.clear();
    clusterIDs.clear();
    boneIDs.clear();
    Mass.clear();
    Vel.clear();
    m_q.clear();
    m_qQ.clear();
    m_Aqq.clear();
    m_AqqQ.clear();
    
}

void ShapeMatch::initShape(const vector<vec> x,
                      const vector<double> mass,
                      const vector<vector<int> > clusIDs,
                      const vector<int> boneID){
    
    if(x.size()!=mass.size()){
        cerr<<"X/mass size difference."<<endl;
        return;
    }reset();
    
    // Append bone cluster
    clusterIDs=clusIDs;
    clusterIDs.push_back(boneID);
    boneIDs = boneID;
    
    Npoint=static_cast<int>(x.size());
    Nclusters=static_cast<int>(clusterIDs.size());
    
    Xorg=x;
    Xnow=Xorg;
    Xnew=Xorg;
    Xtemp=Xorg;
    
    for(int i=0;i<Nclusters;i++){
        vector<vec> goalXC;
        vector<double> massesC;
        for(int j=0;j<clusterIDs[i].size();j++){
            goalXC.push_back(Xnow[clusterIDs[i][j]]);
            massesC.push_back(mass[clusterIDs[i][j]]);
        }
        Xgoal.push_back(goalXC);
        Mass.push_back(massesC);
    }
    vector<vec> vinit(Npoint,zeros<vec>(3));
    Vel=vinit;
    
    m_q.resize(Nclusters);
    m_qQ.resize(Nclusters);
    m_Aqq.resize(Nclusters);
    m_AqqQ.resize(Nclusters);
    third = 1.0/3.0;
    
    free.resize(Npoint);
    for(int i=0;i<Npoint;i++){
        free[i]=true;
    }
    
    if(Npoint<=1) return;
    precompute();
}

void ShapeMatch::relocate(const vector<vec> x){
    
    if(Xorg.size()!=x.size()){
        cerr<<"# new/original points different." <<endl;return;
    }
    Xnow=x;
    Xorg=x;
    Xnew=x;
    Xtemp=x;
    precompute();
}

void ShapeMatch::updateBonePosition(const vector<vec> x){

    for(int i=0;i<x.size();i++){
        Xorg[clusterIDs[Nclusters-1][i]]=x[i];
    }
    //reinitialize final bone cluster
    precompute(Nclusters-1);
}

void ShapeMatch::precompute(){
    
    for(int k=0;k<Nclusters;k++){
        precompute(k);
    }
}

void ShapeMatch::precompute(int k){
    
    vec originalCm=zeros<vec>(3);
    double totalMass=0.;
    for(int i=0;i<Mass[k].size();i++){
        totalMass+=Mass[k][i];
        originalCm+=Xorg[clusterIDs[k][i]]*Mass[k][i];
    }
    originalCm/=totalMass;
    
    mat q=zeros(3,static_cast<int>(clusterIDs[k].size()));
    mat qQ=zeros(9,static_cast<int>(clusterIDs[k].size()));
    mat Aqq=zeros(3,3);
    mat AqqQ=zeros(9,9);
    
    for(int i=0;i<clusterIDs[k].size();i++){
        q.col(i) = Xorg[clusterIDs[k][i]]-originalCm;
        for(int j=0;j<3;j++){
            qQ(j,i)=q(j,i);
            qQ(j+3,i)=q(j,i)*q(j,i);
            qQ(j+6,i)=q(j,i)*q((j+1)%3,i);
        }
        Aqq+=Mass[k][i]*q.col(i)*q.col(i).t();
        AqqQ+=Mass[k][i]*qQ.col(i)*qQ.col(i).t();
    }
    Aqq=pinv(Aqq);
    AqqQ=pinv(AqqQ);
    
    m_q[k]=q;
    m_qQ[k]=qQ;
    m_Aqq[k]=Aqq;
    m_AqqQ[k]=AqqQ;
}

void ShapeMatch::externalForces(){
    
    for(int k=0;k<Nclusters;k++){
        for(int i=0;i<clusterIDs[k].size();i++){
            Xgoal[k][i]=Xorg[clusterIDs[k][i]];
        }
    }
    for(int i=0;i<Npoint;i++){
        Vel[i]+=params.external*params.dt;
        Xnew[i]=Xtemp[i]=Xnow[i]+Vel[i]*params.dt;
    }
}

void ShapeMatch::externalForces(double Fx,double Fy,double Fz){
    vec Fext;
    Fext<<Fx<<Fy<<Fz;
    for(int k=0;k<Nclusters;k++){
        for(int i=0;i<clusterIDs[k].size();i++){
            Xgoal[k][i]=Xorg[clusterIDs[k][i]];
        }
    }
    for(int i=0;i<Npoint;i++){
        Vel[i]+=Fext*params.dt;
        Xnew[i]=Xtemp[i]=Xnow[i]+Vel[i]*params.dt;
    }
}

vector<vector<double> > ShapeMatch::getNormals(void){
    vector<vector<double> > normals(Vel.size());
    for(int i=0;i<normals.size();i++){
        normals[i].resize(3,0.);
    }
    for(int k=0;k<Nclusters;k++){
        vec cm=zeros<vec>(3);
        int h=0.;
        for(int i=0;i<clusterIDs[k].size();i++){
            cm+=Xnew[clusterIDs[k][i]];
            h++;
        }
        cm /= h;
        for(int i=0;i<clusterIDs[k].size();i++){
            vec D = cm-Xnow[clusterIDs[k][i]];
            double d = sqrt(D(0)*D(0)+D(1)*D(1)+D(2)*D(2));
            if(d){
                D/=d;
                normals[clusterIDs[k][i]][0]=D(0);
                normals[clusterIDs[k][i]][1]=D(1);
                normals[clusterIDs[k][i]][2]=D(2);
            }
        }
    }
    return normals;
}


vector<double> ShapeMatch::getCoG(void){
    vec cg=zeros<vec>(3);
    int h=0.;
    for(int i=0;i<Npoint;i++){
        cg+=Xnew[i];
        h++;
    }
    cg /= h;
    vector<double> cog(3,0.);
    cog[0]=cg[0];
    cog[1]=cg[1];
    cog[2]=cg[2];
    return cog;
}

void ShapeMatch::inflate(void){
    
    for(int k=0;k<Nclusters;k++){

        vec cm=zeros<vec>(3);

        int h=0.;
        for(int i=0;i<clusterIDs[k].size();i++){
            cm+=Xnew[clusterIDs[k][i]];
            h++;
        }
        cm /= h;
    
        for(int i=0;i<clusterIDs[k].size();i++){
            vec D = cm-Xnow[clusterIDs[k][i]];
            double d = sqrt(D(0)*D(0)+D(1)*D(1)+D(2)*D(2));
            if(d){
                Vel[clusterIDs[k][i]]+=D/d*params.dt;
                Xnew[clusterIDs[k][i]]=Xtemp[clusterIDs[k][i]]=Xnow[clusterIDs[k][i]]-Vel[clusterIDs[k][i]]*params.dt;
            }
        }
    }
}

/*
void ShapeMatch::inflate(void){
    
    vec cm=zeros<vec>(3);
    int h=0.;
    for(int k=0;k<Nclusters;k++){
        for(int i=0;i<clusterIDs[k].size();i++){
            cm+=Xnew[clusterIDs[k][i]];
            h++;
        }
    }
    cm /= h;

    for(int i=0;i<Npoint;i++){
        vec D = cm-Xnow[i];
        double d = sqrt(D(0)*D(0)+D(1)*D(1)+D(2)*D(2));
        if(d){
            Vel[i]+=D/d*params.dt;
            Xnew[i]=Xtemp[i]=Xnow[i]-Vel[i]*params.dt;
        }
    }
}
 */


void ShapeMatch::projectPositions(){
    
    for(int k=0;k<Nclusters;k++){
        
        vec cm=zeros<vec>(3);   //center of mass
        double totalMass=0.;
        for(int i=0;i<clusterIDs[k].size();i++){
            totalMass+=Mass[k][i];
            cm+=Xnew[clusterIDs[k][i]]*Mass[k][i];
        }
        
        cm /= totalMass;
        mat Apq=zeros(3,3);
        vec p;
        
        for(int i=0;i<clusterIDs[k].size();i++){
            p=Xnew[clusterIDs[k][i]]-cm;
            Apq+=Mass[k][i]*p*m_q[k].col(i).t();
        }
        
        mat R;
        vec S;
        mat u, v;
        svd(u, S, v, Apq);
        R = u * v.t();
        
        switch(params.type){

            case Rigid:{
                
                for(int i=0;i<clusterIDs[k].size();i++){
                    Xgoal[k][i]=R*m_q[k].col(i)+cm;
                    Xnew[clusterIDs[k][i]]+=Xgoal[k][i]-Xnew[clusterIDs[k][i]];
                }
                break;
            }
                
            case Linear:{
                
                mat A=Apq*m_Aqq[k].t();
                if (params.volumeConservation){
                    double det=arma::det(A);
                    if(det!=0.0){
                        det=1.0/pow(fabs(det),third);
                        if(det>3.0) det=3.0;
                        A*=det;
                    }
                }
                mat T=params.beta*A+(1.0-params.beta)*R;
                
                for(int i=0;i<clusterIDs[k].size();i++){
                    /*
                     apply rigid-body dynamics assuming final cluster is bones
                     */
                    if(k==(Nclusters-1)){
                        Xgoal[k][i]=R*m_q[k].col(i)+cm;
                        Xnew[clusterIDs[k][i]]+=Xgoal[k][i]-Xnew[clusterIDs[k][i]];
                    }else{
                        Xgoal[k][i]=T*m_q[k].col(i)+cm;
                        Xnew[clusterIDs[k][i]]+=(Xgoal[k][i]-Xnew[clusterIDs[k][i]])*params.alpha;
                    }
                }
                break;
            }
                
            case Quadratic:{
                
                mat ApqQ=zeros(3,9);
                for(int i=0;i<clusterIDs[k].size();i++){
                    p=Xnew[clusterIDs[k][i]]-cm;
                    ApqQ+=Mass[k][i]*p*m_qQ[k].col(i).t();
                }
                
                mat AQ=ApqQ*m_AqqQ[k];
                mat RQ=join_rows(R,zeros(3,6));
                mat T=params.beta*AQ+(1.0-params.beta)*RQ;
                
                double det=arma::det(T.cols(span(0,2)));

                // This can lead to issues
                if(!params.allowFlip&&det<0.0){
                    T.col(1)=-T.col(1);
                }
                if(params.volumeConservation){
                    if(det!=0.0){
                        det=1.0/pow(fabs(det),third);
                        if(det>3.0) det=3.0;
                        T.cols(span(0,2))*=det;
                    }
                }
                
                for(int i=0;i<clusterIDs[k].size();i++){
                    /*
                     apply rigid-body dynamics assuming final cluster is bones
                     */
                    if(k==(Nclusters-1)){
                        Xgoal[k][i]=R*m_q[k].col(i)+cm;
                        Xnew[clusterIDs[k][i]]+=Xgoal[k][i]-Xnew[clusterIDs[k][i]];
                    }else{
                        Xgoal[k][i]=T*m_qQ[k].col(i)+cm;
                        Xnew[clusterIDs[k][i]]+=(Xgoal[k][i]-Xnew[clusterIDs[k][i]])*params.alpha;
                    }
                }
                break;
            }
            default:
                break;
        }
    }
}

void ShapeMatch::integrate(){
    
    double dt1=1.0/params.dt;
    for(int i=0;i<Npoint;i++){
        if(free[i]){
            Vel[i]=(Xnew[i]-Xnow[i])*dt1;
            Xnow[i]=Xnew[i];
        }
    }
}

void ShapeMatch::iterate(){
    
    externalForces();
    projectPositions();
    integrate();
}