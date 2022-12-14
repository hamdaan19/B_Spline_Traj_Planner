#pragma once

#include <math.h>
#include <nlopt.hpp>
#include <gsl/gsl_poly.h>
#include "B_Spline_Traj_Planner/BSpline.h"

typedef struct {
    Spline3D H;
} constraint_data; 

class GFunc {
    /* class entity for G function which is a direct relationship 
    between my time 't' domain and tau domain */

    public:
        double alpha1, alpha2, alpha3, alpha4, alpha5;
        double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4,
               e1, e2, e3, e4, f1, f2, f3, f4, g1, g2, g3, g4, h1, h2, h3, h4; 

        int degree; // Degree of the polynomial. For cubic, degree = 3
        double totalSegs;
        std::vector<double> segIntervals; 

        std::vector<double> alpha;

        GFunc(std::vector<double> alpha, int degree);

        double G(double tau);
        double invertedG(double t);
        double gradG(double tau);
        double grad2G(double tau);
        double grad3G(double tau);
        double minimumGradG();

        double maximumPositionXGradient(Spline3D H);
        double maximumPositionYGradient(Spline3D H);
        double maximumPositionZGradient(Spline3D H);

        double maximumVelocityXGradient(Spline3D H);
        double maximumVelocityYGradient(Spline3D H);
        double maximumVelocityZGradient(Spline3D H);

        double sumOfJerks(Spline3D H);


    private:
        static constexpr double dtau = 1e-2; 
        int checkWhichSegment(double tau); 
        int checkWhichTimeSegment(double t);
};

