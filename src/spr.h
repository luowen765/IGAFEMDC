// Copyright (c) 2024.
// This file is part of the IGAFEMDC program. IGAFEMDC is free software, you can
// redistribute it and/or modify it under the terms of the BSD-3 license. See
// file LICENSE for details.
// For more information and source code availability, please visit
// https://github.com/luowen765/IGAFEMDC.
// @Author: Lewen liu; Jingtan Tang, Zhengguang liu, Zhengyong Ren.

#ifndef _Point_ren_H
#define _Point_ren_H

// C++ includes
#include <cassert>
#include <cmath>
#include <iostream>
// Local includes

class Point_ren {
public:
  Point_ren(const double x = 0., const double y = 0., const double z = 0.);
  Point_ren(const Point_ren &p);
  virtual ~Point_ren();

  friend std::ostream &operator<<(std::ostream &os, const Point_ren &Point_ren);

  // Overload operators.
  Point_ren &operator=(const Point_ren &p);
  double operator()(const unsigned int i) const;
  double &operator()(const unsigned int i);
  Point_ren operator+(const Point_ren &v) const;
  Point_ren operator-(const Point_ren &v) const;
  Point_ren operator*(const double a) const;
  double operator*(const Point_ren &v) const;
  bool operator==(const Point_ren &v) const;
  bool operator<(const Point_ren &v) const;
  void operator=(const double a);

  // Math functions
  Point_ren cross(const Point_ren &v) const;
  Point_ren unit() const;
  double size() const;
  void zero();
  double *get_xyz() { return _coords; }

protected:
  double _coords[3];
};

//-----------------------------------------------------------------------
inline Point_ren::Point_ren(const double x, const double y, const double z) {
  _coords[0] = x;
  _coords[1] = y;
  _coords[2] = z;
}

inline Point_ren::Point_ren(const Point_ren &p) {
  for (unsigned int i = 0; i < 3; i++)
    _coords[i] = p._coords[i];
}

inline Point_ren &Point_ren::operator=(const Point_ren &p) {
  _coords[0] = p._coords[0];
  _coords[1] = p._coords[1];
  _coords[2] = p._coords[2];

  return (*this);
}

inline Point_ren::~Point_ren() {
  // no space is allocated by the new operator.
  // so,there is no delete [].
}

inline double Point_ren::operator()(const unsigned int i) const {
  assert(i < 3);
  return _coords[i];
}

inline double &Point_ren::operator()(const unsigned int i) {
  assert(i < 3);
  return _coords[i];
}

inline Point_ren Point_ren::operator+(const Point_ren &p) const {
  return Point_ren(_coords[0] + p._coords[0], _coords[1] + p._coords[1],
                   _coords[2] + p._coords[2]);
}

inline Point_ren Point_ren::operator-(const Point_ren &p) const {
  return Point_ren(_coords[0] - p._coords[0], _coords[1] - p._coords[1],
                   _coords[2] - p._coords[2]);
}

inline Point_ren Point_ren::operator*(const double factor) const {
  return Point_ren(_coords[0] * factor, _coords[1] * factor,
                   _coords[2] * factor);
}

inline double Point_ren::operator*(const Point_ren &p) const {
  return (_coords[0] * p(0) + _coords[1] * p(1) + _coords[2] * p(2));
}

inline bool Point_ren::operator==(const Point_ren &rhs) const {
  return ((fabs(_coords[0] - rhs._coords[0]) +
           fabs(_coords[1] - rhs._coords[1]) +
           fabs(_coords[2] - rhs._coords[2])) < 3 * 1e-6);
}

inline bool Point_ren::operator<(const Point_ren &rhs) const {
  // First we assume (this)<rhs true
  if (*this == rhs)
    return false;
  if ((*this)(0) < rhs(0))
    return true; //  <
  else if ((*this)(0) > rhs(0))
    return false;                                  //  >
  else if (std::abs((*this)(0) - rhs(0)) < 1e-6) { // vx=rhsx
    if ((*this)(1) < rhs(1))
      return true;
    else if ((*this)(1) > rhs(1))
      return false;
    else if (std::abs((*this)(1) - rhs(1)) < 1e-6) { // vy=rhsy
      if ((*this)(2) < rhs(2))
        return true;
      else if ((*this)(2) > rhs(2))
        return false;
    }
  }
  return false;
}

inline void Point_ren::operator=(const double a) {
  _coords[0] = a;
  _coords[1] = a;
  _coords[2] = a;
}

inline Point_ren Point_ren::cross(const Point_ren &p) const {
  return Point_ren(_coords[1] * p._coords[2] - _coords[2] * p._coords[1],
                   -_coords[0] * p._coords[2] + _coords[2] * p._coords[0],
                   _coords[0] * p._coords[1] - _coords[1] * p._coords[0]);
}

inline Point_ren Point_ren::unit() const {
  const double length = size();
  return Point_ren(_coords[0] / length, _coords[1] / length,
                   _coords[2] / length);
}

inline double Point_ren::size() const {
  double value = std::pow(_coords[0], 2) + std::pow(_coords[1], 2) +
                 std::pow(_coords[2], 2);

  return std::sqrt(value);
}

inline void Point_ren::zero() {
  _coords[0] = 0.;
  _coords[1] = 0.;
  _coords[2] = 0.;
}

inline std::ostream &operator<<(std::ostream &os, const Point_ren &Point_ren) {
  os << Point_ren(0) << "\t" << Point_ren(1) << "\t" << Point_ren(2);
  return os;
}

#endif // _Point_ren__H

#ifndef _SPR_H
#define _SPR_H

// C++ include
#include <set>
#include <string>
#include <vector>

#include "dc.h"
#include "mfem.hpp"
#include "parahandler.h"
//---------------------------------------------------------------------------
// Class SPR
class SPR {
public:
  // Destructor.
  ~SPR(){};

  SPR(ParaHandler *para_handler_, ParMesh *pmesh, ParFiniteElementSpace *pfes,
      const unsigned int s = 10); // size of patch.

  void gradient_patch_recovery(std::vector<DenseVector> &Ui, DenseVector &W,
                               std::vector<DenseMatrix> &cond_att,
                               Vector &local_err);

protected:
  // Get the element patch "patch" surrounding the element "elem".
  void get_element_patch(unsigned int elem,                // a element
                         std::set<unsigned int> &patch_id, // element patch.
                         const unsigned int s,             // size of patch.
                         //  std::vector<Matrix3D> &sigma
                         std::vector<DenseMatrix> &cond_att);
  // Get the polynomial terms with the same order p of finite-element
  // "elem" in the patch surrounding element "elem". Its value is adopted at
  // point "p".
  void get_patch_polynomial_c(std::vector<double> &c, // polynomial terms
                              const Point_ren &p);    // sampling points.
  // Get the polynomial parameters matrix for 3-component in the element
  // patch at all sampling points. Herein, a least-square fit method will be
  // adopted.
  void get_patch_b(std::set<unsigned int> &patch_id,
                   std::vector<DenseVector> &solution,
                   std::vector<Matrix43D> &b);
  // Calculate Pi_{k} in equation (23)
  void compute_PI_k(std::set<unsigned int> &patch_id,
                    std::vector<DenseVector> &solution,
                    std::vector<Matrix43D> &a,
                    // std::vector<Matrix3D> &sigma,
                    std::vector<DenseMatrix> &cond_att,
                    std::map<unsigned int, double> &energy_norm);
  // compute the C 4x4 symmetrical matrix
  void compute_C(std::set<unsigned int> &patch_id, Matrix4D &C);
  void solve_a(Matrix4D &C, std::vector<Matrix43D> &b,
               std::vector<Matrix43D> &a);
  void compute_P_value(Point_ren &r, Matrix43D &a, Matrix13D &P);

private:
  ParMesh *_pmesh;
  ParFiniteElementSpace *_pfes;
  ParaHandler *_para_handler;
  // The size of element patch.
  unsigned int _s;
};

#endif // SPR_H
