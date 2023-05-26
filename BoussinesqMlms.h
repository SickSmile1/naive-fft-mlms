/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef BOUSSINESQMLMS_H_
#define BOUSSINESQMLMS_H_

#include "Boussinesq.h"

matrix initializeStylusArray(int t); // NOLINT

void initializeStack(matrix &st, const int t, const matrix Ip, // NOLINT
                    const matrix kM, int grid, // NOLINT
                    std::vector<matrix>& pfVec,
                    std::vector<matrix>& cDVec);

// most outer loop, for we need a n*n loop for the algorithm
void naiveCalculation(matrix &Ic, const matrix &Pa, double cell_size); // NOLINT

// inner loop calling the calulation for every n
inline double calcBoussinesq(double a, double b, double x, double y);

inline double calcBoussinesq(int i, int j, double dxc, double dyc,
    double dxf, double dyf);

// actual calculation
double calc_displacement(const matrix &pressure,
                const matrix &Ic,
                double y, double x,
                double a, double b,
                double cell);

void calc_displacement(const matrix &pF, double cS, double fS,
                         matrix &cD);  // NOLINT

void calcCoarsePressure(std::vector<matrix> &pFVec, // NOLINT
                        const matrix& st);

// coarse pressure array calculation
// void calc_coarse_pressure(const matrix &fP, const matrix &st,
//                           matrix &cP, std::size_t t); // NOLINT


// calculate the sizes for the coarse grids
void prepareCoarseSizes(std::vector<std::size_t> &gridLen1, // NOLINT
                std::vector<std::size_t> &gridLen2, // NOLINT
                std::vector<std::size_t> &coarseGrid, // NOLINT
                const std::size_t shape,
                const std::size_t t);


bool boundaryCheck(matrix &m, int i, int j); // NOLINT
bool boundaryCheck(const matrix &m, int i, int j); // NOLINT

void correctionSteps(matrix& cC, const matrix& st, int mc, // NOLINT
                    int t, double fineSize, double halfSize);


void applyCorrection(matrix &coarseDisplacement, const matrix cC, // NOLINT
                     const matrix Ip, int t);

double correctionHelper(const matrix& cC, const matrix& Ip, int t, // NOLINT
                        int i, int j);

void interpolateGrid(matrix &nextDisplacement, // NOLINT
                     const matrix coarseDisplacement, // NOLINT
                     const matrix st);

void secondCorrectionStep(const matrix& st, // NOLINT
                          double hS, const matrix& pF, // NOLINT
                          matrix& cD, const std::vector<matrix> &cCVec);  // NOLINT


void createCorrectionArrays(std::vector<matrix> &cCVec, // NOLINT
                            const matrix &st, double hS, // NOLINT
                            double fineSize);

#endif  // BOUSSINESQMLMS_H_

