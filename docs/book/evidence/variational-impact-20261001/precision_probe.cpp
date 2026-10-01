// Independent probes of event-location convergence and orbital momentum in oblique impacts.
#include "rigid/RigidWorld.h"
#include "core/Format.h"
#include <iostream>
#include <cmath>
using namespace rf;
void setup(RigidWorld& w) {
    w.params.sleeping=w.params.collideWithDomain=false;
    w.params.linearDamping=w.params.angularDamping=w.params.rollingResistance=0;
    w.params.iterations=40;
}
int add(RigidWorld& w, Vector3 p, float mass) {
    int id=w.addSphere(p,0.1f,1,Vector3(1));auto& b=w.bodies()[id];
    b.mass=mass;b.invMass=1/mass;b.invInertiaLocal=Vector3(250/mass);b.updateInertia();
    b.restitution=1;b.friction=b.staticFriction=0;return id;
}
double energy(const RigidWorld& w) {
    double e=0;for(const auto& b:w.bodies())if(b.invMass>0)for(int k=0;k<3;++k)
        e+=b.mass*(0.5*double(b.vel[k])*b.vel[k]-double(w.params.gravity[k])*b.pos[k]);
    return e;
}
int convergence() {
    for(float tol:{0.002f,0.001f,0.0005f,0.00025f}) {
        RigidWorld w;setup(w);w.params.ccdTolerance=tol;
        w.addBox({0,-5,0},{25,5,25},Quaternion(),0,Vector3(1));int id=add(w,{0,5,0},1);
        const double g=-w.params.gravity.y,hit=std::sqrt(2*4.9/g),e0=energy(w);
        double t=0,qsq=0,worst=0;size_t events=0;
        for(int i=0;i<2000;++i) {
            if(!w.tryVariationalStep(0.01f))return 2;
            t+=0.01f;events+=w.lastStepResult().impactEvents;
            double tau=std::fmod(t-hit,2*hit);
            double exact=t<hit?5-0.5*g*t*t:0.1+g*hit*tau-0.5*g*tau*tau;
            qsq+=std::pow(w.bodies()[id].pos.y-exact,2);worst=std::max(worst,std::fabs(energy(w)/e0-1));
        }
        std::cout<<"{\"tolerance_m\":"<<numberText(tol)<<",\"position_rmse_m\":"<<numberText(float(std::sqrt(qsq/2000)))
                 <<",\"max_relative_energy_error\":"<<numberText(float(worst))<<",\"events\":"<<events<<"}\n";
    }
    return 0;
}
int oblique() {
    float maxL=0,maxP=0;double maxE=0;size_t events=0;
    for(int trial=0;trial<40;++trial) {
        RigidWorld w;setup(w);w.params.gravity=Vector3(0);
        const float y=0.005f+trial*0.002f;
        int a=add(w,{-0.5f,y,0},1),b=add(w,{0.5f,-y,0},2);
        w.bodies()[a].vel={2,0,0};w.bodies()[b].vel={-1,0,0};
        const auto momentum=[&](){return w.bodies()[a].vel+2*w.bodies()[b].vel;};
        const auto angular=[&](){return cross(w.bodies()[a].pos,w.bodies()[a].vel)+2*cross(w.bodies()[b].pos,w.bodies()[b].vel);};
        const Vector3 p=momentum(),l=angular();const double e=energy(w);
        for(int k=0;k<50;++k){if(!w.tryVariationalStep(0.01f))return 3;events+=w.lastStepResult().impactEvents;}
        maxP=std::max(maxP,length(momentum()-p));maxL=std::max(maxL,length(angular()-l));
        maxE=std::max(maxE,std::fabs(energy(w)/e-1));
    }
    std::cout<<"{\"oblique_cases\":40,\"events\":"<<events<<",\"max_momentum_error\":"<<numberText(maxP)
             <<",\"max_orbital_angular_momentum_error\":"<<numberText(maxL)
             <<",\"max_relative_energy_error\":"<<numberText(float(maxE))<<"}\n";
    return events!=40||maxP>1e-5f||maxL>1e-5f||maxE>1e-5;
}
int main(){const int code=convergence();return code?code:oblique();}
