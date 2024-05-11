
// Copyright (c) 2024.
// This file is part of the IGAFEMDC program. IGAFEMDC is free software, you can
// redistribute it and/or modify it under the terms of the BSD-3 license. See
// file LICENSE for details.
// For more information and source code availability, please visit
// https://github.com/luowen765/IGAFEMDC.
// @Author: Lewen liu; Jingtan Tang, Zhengguang liu.

/*
 * @Description:
  This class does post-processing for the solutions of FEM. The total potential
can be obtained by adding the primary singular potential. Then apparent
resistivity can be computed on different configrations. 
 */

#ifndef _POST_H
#define _POST_H

#include "dc.h"
#include "mfem.hpp"
#include "mfemplus.h"
#include "parahandler.h"
#include <fstream>
#include <string>
#include <vector>

using namespace mfem;

class Post {
public:
  ParMesh *pmesh;
  ParFiniteElementSpace *pfes;
  std::vector<ParGridFunction *> up;
  ParaHandler *para_handler;
  std::vector<Array<int>> local_sites_tets;
  std::vector<Array<double>> local_Up;
  std::vector<Array<double>> Ex;
  std::vector<Array<double>> Ey;
  std::vector<Array<double>> Ez;
  std::vector<Array<double>> Jx;
  std::vector<Array<double>> Jy;
  std::vector<Array<double>> Jz;
  std::vector<Array<double>> Ut;

public:
  Post(ParMesh *pmesh_, ParFiniteElementSpace *pfes_,
       std::vector<ParGridFunction *> &up_, ParaHandler *para_handler_);

  ~Post();

public:
  // Post-processing and save the solution
  void main_post_process(ParaHandler *para_handler, std::string amr,
                         std::vector<DenseMatrix> &sigma0,
                         std::vector<DenseMatrix> &cond_att);

  // Post-processing
  void post_process(std::vector<DenseMatrix> &cond_att,
                    std::vector<DenseMatrix> &sigma0);
  // Save results of global sites in rank 0 (root rank)
  void main_save_local_mesh(std::string amr_, Meshplus *pmesh,
                            Vector &local_eta,
                            std::vector<DenseMatrix> &att_cond);
  // compute total potential Ut(c,m) and apparent resistivity rho, print them
  void read_compute_survery_mode(std::string survery_mode_file,
                                 std::string output_file);

  void compute_oneSource_U_E_J_averagely(std::vector<DenseMatrix> &cond_att,
                                         std::vector<DenseMatrix> &sigma0);

  void saveGridFunction(std::string amr_, std::string fn, Meshplus *pmesh,
                        ParGridFunction *field);
};

#endif // _POST_H
