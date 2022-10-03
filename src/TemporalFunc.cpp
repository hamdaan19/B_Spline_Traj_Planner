#include <iostream>
#include <nlopt.hpp>
#include <cmath>

#include "B_Spline_Traj_Planner/TemporalFunc.h"

GFunc::GFunc(std::vector<double> alpha, int degree) {

    this->alpha = alpha; 
    this->degree = degree; 

    totalSegs = alpha.size()/(degree+1);

    for (double interval = 0; interval <= 1; interval += (1/totalSegs)){
        segIntervals.push_back(interval);
    }

}

int GFunc::checkWhichSegment(double tau) {
    int seg; 
    for (int idx = 1; idx < segIntervals.size(); idx++){
        if (tau <= segIntervals.at(idx)){
            seg = idx;
            break;  
        }
    }
    return seg; 
}

int GFunc::checkWhichTimeSegment(double t){
    int seg;
    for (int idx = 1; idx < segIntervals.size(); idx++){
        if (t <= G(segIntervals.at(idx))){
            seg = idx;
            break;
        }
    }
    return seg;
}

double GFunc::invertedG(double t){
    int seg = checkWhichTimeSegment(t);
    int k = (degree+1)*(seg-1);

    double coeff[degree+1]; 
    double roots[degree]; 

    for (int i = 0; i < degree; i++){
        coeff[i] = alpha[k+i]/alpha[k];
    }
    coeff[3] = (alpha[k+3]-t)/alpha[k]; 

    gsl_poly_solve_cubic(coeff[1], coeff[2], coeff[3], &roots[0], &roots[1], &roots[2]);
    
    return roots[0]; 
    
}

double GFunc::G(double tau){ //tau exists between 0 and 1

    double t = 0; 
    int seg = checkWhichSegment(tau);

    int i = (degree+1)*(seg-1);
    for (int k = 0; k <= degree; k++){
        t += alpha[i+k]*pow(tau, degree-k);
    }

    return t;
}

double GFunc::gradG(double tau){
    double t = 0; 
    int seg = checkWhichSegment(tau);

    int i = (degree+1)*(seg-1);
    for (int k = 0; k <= degree-1; k++){
        t += (degree-k)*alpha[i+k]*pow(tau, degree-k-1);
    }

    return t; 
}

double GFunc::grad2G(double tau){
    double t = 0;
    int seg = checkWhichSegment(tau);

    int i = (degree+1)*(seg-1);
    for (int k = 0; k<=degree-2; k++){
        t += (degree-k-1)*(degree-k)*alpha[i+k]*pow(tau, degree-k-2);  
    }

    return t; 
}

double GFunc::grad3G(double tau){
    double t = 0; 
    int seg = checkWhichSegment(tau);
    
    int i = (degree+1)*(seg-1);
    for (int k = 0; k<=degree-3; k++){
        t += (degree-k-2)*(degree-k-1)*(degree-k)*alpha[i+k]*pow(tau, degree-k-3);
    }

    return t; 
}

double GFunc::minimumGradG(){ 
    // Incrementally checks the value of dG with step size of dtau.
    double minGrad = gradG(0);
    for (double tau = 0; tau <= 1; tau += dtau){
        if (gradG(tau) < minGrad) {
            minGrad = gradG(tau);
        }
    }
    return minGrad; 
}

double GFunc::maximumPositionXGradient(Spline3D H){

    double maxVel = 0;
    
    for (double tau = 0; tau <= 1; tau += dtau){
        auto deriv = H.derivatives(tau, 1);

        double velX = deriv(0, 1) * (1/gradG(tau));

        if (abs(velX) > maxVel){
            maxVel = abs(velX);
        }
    }

    return maxVel;
}

double GFunc::maximumPositionYGradient(Spline3D H){

    double maxVel = 0;
    
    for (double tau = 0; tau <= 1; tau += dtau){
        auto deriv = H.derivatives(tau, 1);

        double velY = deriv(1, 1) * (1/gradG(tau));

        if (abs(velY) > maxVel){
            maxVel = abs(velY);
        }
    }

    return maxVel;
}

double GFunc::maximumPositionZGradient(Spline3D H){

    double maxVel = 0;
    
    for (double tau = 0; tau <= 1; tau += dtau){
        auto deriv = H.derivatives(tau, 1);

        double velZ = deriv(2, 1) * (1/gradG(tau));

        if (abs(velZ) > maxVel){
            maxVel = abs(velZ);
        }
    }

    return maxVel;
}

double GFunc::maximumVelocityXGradient(Spline3D H){
    double maxAcc = 0;

    for (double tau = 0; tau <= 1 - 2*dtau; tau += dtau){
        auto deriv = H.derivatives(tau, 2);

        double accX = ( 1/gradG(tau) ) * ( (deriv(0, 2)/gradG(tau)) - ( (deriv(0, 1)*grad2G(tau))/pow(gradG(tau), 2) ) );
        
        if (abs(accX) > maxAcc) {
            maxAcc = abs(accX);
        } 
    }

    return maxAcc;
}

double GFunc::maximumVelocityYGradient(Spline3D H){
    double maxAcc = 0;

    for (double tau = 0; tau <= 1 - 2*dtau; tau += dtau){
        auto deriv = H.derivatives(tau, 2);

        double accY = ( 1/gradG(tau) ) * ( (deriv(1, 2)/gradG(tau)) - ( (deriv(1, 1)*grad2G(tau))/pow(gradG(tau), 2) ) );

        if (abs(accY) > maxAcc) {
            maxAcc = abs(accY);
        } 
    }

    return maxAcc;
}

double GFunc::maximumVelocityZGradient(Spline3D H){
    double maxAcc = 0;

    for (double tau = 0; tau <= 1 - 2*dtau; tau += dtau){
        auto deriv = H.derivatives(tau, 2);

        double accZ = ( 1/gradG(tau) ) * ( (deriv(2, 2)/gradG(tau)) - ( (deriv(2, 1)*grad2G(tau))/pow(gradG(tau), 2) ) );

        if (abs(accZ) > maxAcc) {
            maxAcc = abs(accZ);
        } 
    }

    return maxAcc;
}



double GFunc::sumOfJerks(Spline3D H){
    double sum = 0;

    for (double tau = 0; tau <= 1; tau += dtau){
        auto deriv = H.derivatives(tau, 3);

        auto H = deriv.col(0);
        auto d1H = deriv.col(1);
        auto d2H = deriv.col(2);
        auto d3H = deriv.col(3);

        double D = pow(gradG(tau), 3); 
        double gradD = 3*pow(gradG(tau), 2)*grad2G(tau);
        Eigen::Vector3d N = gradG(tau)*d2H - d1H*grad2G(tau);
        Eigen::Vector3d gradN = gradG(tau)*d3H - d1H*grad3G(tau);

        Eigen::Vector3d jerk = (1/gradG(tau)) * (D*gradN - N*gradD)/pow(D, 2); 

        sum += abs(jerk.squaredNorm());
 
    }

    return sum;
}




