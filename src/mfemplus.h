// Copyright (c) 2024.
// This file is part of the IGAFEMDC program. IGAFEMDC is free software, you can
// redistribute it and/or modify it under the terms of the BSD-3 license. See
// file LICENSE for details.
// For more information and source code availability, please visit
// https://github.com/luowen765/IGAFEMDC.
// @Author: Lewen liu; Jingtan Tang, Zhengguang liu.


// This class inherits MFEM's class in mesh/mesh.hpp. I do this for extending
// the print function This file is part of the MFEM library. For more
// information and source code availability visit https://mfem.org. MFEM is free
// software; you can redistribute it and/or modify it under the terms of the
// BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "mfem.hpp"
using namespace mfem;

#ifndef _MESHPLUS_H
#define _MESHPLUS_H

// #include "mfem.hpp"
// using namespace mfem;

class Meshplus : public Mesh {

  // class Meshplus : public ParMesh {

public:
  Meshplus();
  Meshplus(ParMesh *pmesh_);
  Meshplus(Mesh *mesh_);
  void PrintVTK_eta(std::ostream &out, Vector &other_attri);
  void write_out_node_info_vtk(std::ostream &out, Vector &node_data);
};

#endif // _MESHPLUS_H





/*
 * @Description:
The following classes are the direct copies and modifies of MFEM's classes
in fem/coefficient.hpp. I just did a minor modifications to compute different
coeffients for assembling equations. For more information and source code
availability, please visit https://github.com/luowen765/3DDCAF.
 */

#ifndef MFEM_COEFFICIENTS_H
#define MFEM_COEFFICIENTS_H

#include "dc.h"
#include "mfem.hpp"
using namespace mfem;

class LVCoefficient : public MatrixCoefficient {
private:
  std::vector<DenseMatrix> mat;
  int temp_count;

public:
  /// Construct a coefficeint using matrix @a m for left side volume integral
  LVCoefficient(std::vector<DenseMatrix> &m)
      : MatrixCoefficient(m[0].Height(), m[0].Width()), mat(m) {}
  using MatrixCoefficient::Eval;
  /// Evaluate the matrix coefficient at @a ip.
  virtual void Eval(DenseMatrix &M, ElementTransformation &T,
                    const IntegrationPoint &ip);
};

class LGamma1Coefficient : public Coefficient {
private:
  std::vector<DenseMatrix> mat;
  Vector x_c;
  ParMesh *pmesh;
  bool boundaryFlag = true;

public:
  /// Construct a double coefficient for the surface integral in distant
  /// boundary face
  LGamma1Coefficient(std::vector<DenseMatrix> &m, Vector sources_center,
                     ParMesh *pmesh_)
      : mat(m), x_c(sources_center), pmesh(pmesh_) {}
  using Coefficient::Eval;
  /// Evaluate the coefficient at @a ip.
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};

class RVCoefficient : public VectorCoefficient {
private:
  std::vector<DenseMatrix> mat;
  DenseMatrix sigma0;
  Vector x_i;
  // ParMesh *pmesh;

public:
  /// Construct a coeffient  using vector @a V for right side volume integral
  RVCoefficient(std::vector<DenseMatrix> &m, DenseMatrix sigma0_,
                Vector source_i)
      : VectorCoefficient(3), mat(m), sigma0(sigma0_), x_i(source_i)
  // ,pmesh(pmesh_)
  {}
  using VectorCoefficient::Eval;
  /// Evaluate the vector coefficient at @a ip.
  virtual void Eval(Vector &V, ElementTransformation &T,
                    const IntegrationPoint &ip);
};

class RGamma0Coefficient : public Coefficient {
private:
  DenseMatrix sigma0;
  Vector x_i;
  ParMesh *pmesh;
  bool boundaryFlag = true;

public:
  /// Construct a coeffient for right side volume integral
  RGamma0Coefficient(DenseMatrix sigma0_, Vector source_i, ParMesh *pmesh_)
      : sigma0(sigma0_), x_i(source_i), pmesh(pmesh_) {}
  using Coefficient::Eval;
  /// Evaluate the double coefficient at @a ip.
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};

class RGamma1Coefficient : public Coefficient {
private:
  std::vector<DenseMatrix> mat;
  DenseMatrix sigma0;
  Vector x_i;
  ParMesh *pmesh;
  bool boundaryFlag = true;

public:
  /// Construct a coeffient for right side volume integral
  RGamma1Coefficient(std::vector<DenseMatrix> &m, DenseMatrix sigma0_,
                     Vector source_i, ParMesh *pmesh_)
      : mat(m), sigma0(sigma0_), x_i(source_i), pmesh(pmesh_) {}
  using Coefficient::Eval;
  /// Evaluate the double coefficient at @a ip.
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};

#endif