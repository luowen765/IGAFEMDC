

// Copyright (c) 2024.
// This file is part of the IGAFEMDC program. IGAFEMDC is free software, you can
// redistribute it and/or modify it under the terms of the BSD-3 license. See
// file LICENSE for details.
// For more information and source code availability, please visit
// https://github.com/luowen765/IGAFEMDC.
// @Author: Lewen liu; Jingtan Tang, Zhengguang liu.

#ifndef _EM_H
#define _EM_H

#include <Eigen/Dense>  // Linear Algebra Lib.
#include <Eigen/Sparse> // Sparse Lib
#include <Eigen/StdVector>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>

#include "mfem.hpp"
#include <string>
#include <sys/resource.h>
#include <unistd.h>
using namespace mfem;

namespace DC {
typedef Eigen::VectorXd DenseVector;          // double
typedef Eigen::Matrix<double, 3, 3> Matrix3D; // 3*3
typedef Eigen::Matrix<double, 4, 3> Matrix43D;
typedef Eigen::Matrix<double, 1, 3> Matrix13D;
typedef Eigen::Matrix<double, 1, 4> Matrix14D;
typedef Eigen::Matrix<double, 4, 4> Matrix4D;

static const double PI = 3.1415926535897932384626433832795;
// Tolerance used in our package.
static const double TOLERANCE = 1e-6;
} // namespace DC
using namespace DC;

// Calculate the length between two points: coord1 and coord2
double length_two_point(double *coord1, double *coord2);
double length_two_point(Vector &coord1, Vector &coord2);

void remove_duplicate_vertex(std::vector<mfem::Vertex> &all_point);

double B_value(Vector p, Vector pi, DenseMatrix sigma);

double Beta_value(Vector p, Vector pi, DenseMatrix sigma, Vector normal);

double U_i_s(Vector p, Vector pi, DenseMatrix sigma0);

void gradient_U_i_s(mfem::Vector p, mfem::Vector pi, mfem::DenseMatrix sigma0,
                    mfem::Vector &deltaUs);

void compute_normal(ElementTransformation &T, const IntegrationPoint &ip,
                    Vector &normal);

void get_tetId_bdr(ParMesh *pmesh, int ElementNo, int &tet_id);

// whether a==b, note width == height for a and b
bool equalDM(DenseMatrix a, DenseMatrix b);
bool equalVector(Vector a, Vector b);

void Mult_DenseMatrix3(DenseMatrix &a, DenseMatrix &b, DenseMatrix &c);

#endif // _EM_H
