

// Copyright (c) 2024.
// This file is part of the IGAFEMDC program. IGAFEMDC is free software, you can
// redistribute it and/or modify it under the terms of the BSD-3 license. See
// file LICENSE for details.
// For more information and source code availability, please visit
// https://github.com/luowen765/IGAFEMDC.
// @Author: Lewen liu; Jingtan Tang, Zhengguang liu.

#ifndef _GOAFEM_H
#define _GOAFEM_H

#include "dc.h"
#include "mfem.hpp"
#include "parahandler.h"

#include <mpi.h>
#include <vector>

using namespace std;
using namespace mfem;

class GOAFEM {
public:

  std::vector<DenseMatrix> cond_att;
  std::vector<DenseMatrix> sigma0;

  ParaHandler *para_handler;

  ParMesh *pmesh;

  FiniteElementCollection *pfec;

  ParFiniteElementSpace *pfes;
  BilinearFormIntegrator *integ;

  ParBilinearForm *a;
  // Primal problem a(u,v) = f(v)
  std::vector<ParGridFunction *> up;
  // Dual problem a(w,v) = l(v)
  ParGridFunction *w;
  std::vector<ParLinearForm *> f;
  // dual problem linear form
  ParLinearForm *l;

  Array<int> gamma0_bdr;
  Array<int> gamma1_bdr;

  MatrixCoefficient *LVolume_coef;
  Coefficient *LGamma1_coef;
  std::vector<VectorCoefficient *> RVolume_coef;
  std::vector<Coefficient *> RGamma0_coef;
  std::vector<Coefficient *> RGamma1_coef;
  OperatorHandle A;
  std::vector<Vector> U, F;
  Vector W, L;
  Vector local_err;
  Array<int> local_sources_tets;

  Array<int> global_initial_sources_tets;

  Array<int> local_dual_tets;

  int order;

public:
  GOAFEM(std::vector<DenseMatrix> &att_cond, ParaHandler &para_handler_,
         ParMesh *pmesh_);
  ~GOAFEM();

public:
  void initialize();

  /*   Setup coefficient for left volume integral
    Setup coefficient for left distant surface integral
    Setup coefficient for right volume integral
    Setup coefficient for right distant surface integral
    Setup coefficient for air-earth surface integral */
  void setup_integral_coefficient();
  void set_bdr();

  // Problem size
  HYPRE_Int get_problem_dofs();
  void print_problem_size();
  void assemble_dual_linearform();
  void solve_with_pcg(OperatorHandle &A, Vector &X, Vector &B, std::string pre,
                      std::string problem_type);

  // solve Boundary Value Problem with different solvers
  void solve(OperatorHandle &A, std::vector<Vector> &X, std::vector<Vector> &B,
             std::string problem_type);

  void solve_primal_problem();

  void solve_dual_problem();

  // Update all objects based on pmesh
  void set_sigma0();
  // Estimating error
  void error_estimating();

  void refine_mesh();
};
#endif // _GOAFEM_H
