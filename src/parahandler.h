/*
 * @Description: 
 * @version: 
 * @Author: luowen
 * @Date: 2024-03-24 22:16:38
 * @LastEditTime: 2024-05-06 21:44:20
 * @FilePath: /IGAFEMDC/src/parahandler.h
 */
// Copyright (c) 2024.
// This file is part of the IGAFEMDC program. IGAFEMDC is free software, you can
// redistribute it and/or modify it under the terms of the BSD-3 license. See
// file LICENSE for details.
// For more information and source code availability, please visit
// https://github.com/luowen765/IGAFEMDC.
// @Author: Lewen liu; Jingtan Tang, Zhengguang liu.

/*
 * @Description:
 A class for managing input information.
 */

#ifndef _PARAHANDLER_H
#define _PARAHANDLER_H

#include "dc.h"
#include "mfem.hpp"
#include <fstream>
#include <map>
#include <math.h>
#include <string>
#include <vector>

using namespace DC;
using namespace mfem;
class ParaHandler {
public:
  int source_number;
  std::vector<Vertex> sources;
  std::vector<std::vector<Vertex>> sites;
  std::vector<Vertex> s_plus_m;
  std::string survey_input;

  std::string model_parameters_file;
  std::string marker_type;
  int n_regions;
  std::vector<int> marker_vec;
  std::map<int, DenseMatrix> region2conductivity;

  int maxit;         // max AMR iterations
  long int max_dofs; // max dofs
  int amr_type;            // the method of choosing the refined mesh
  double beta;       // marking parameter
  std::string linear_solver;

  int save_amr_mesh;
  std::string mesh_file;

  std::string if_compute_P2;
  std::string if_compute_J;
  std::string if_use_improved_AMR;


public:
  ParaHandler(char *model_file);
  ~ParaHandler();

public:
  void skip_comments(std::istream &in, std::vector<std::string> &para_str_vec,
                     std::string comment_str = "#");
  // read model information from file
  void read_model_info(char *model_file);
  // calculate the tensor conductivity
  DenseMatrix cal_conductivity(Vector &main_cond, Vector &_angle);
  // get element conductivity according to marker (attribute)
  DenseMatrix get_elem_conductivity(int marker);
  // get center coordinates of sources
  Vector get_sources_center();

  void preProcess();
  void find_point_tets(ParMesh *pmesh, std::vector<Vertex> &points,
                       Array<int> &find_tets, std::string find_points_by);
};
#endif // _PARAHANDLER_H
