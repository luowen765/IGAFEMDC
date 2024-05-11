// Copyright (c) 2023.
// This file is part of the IGAFEMDC program. IGAFEMDC is free software, you can
// redistribute it and/or modify it under the terms of the BSD-3 license. See
// file LICENSE for details.
// For more information and source code availability, please visit
// https://github.com/luowen765/IGAFEMDC.
// @Author: Lewen liu; Jingtan Tang, Zhengguang liu.

#include "dc.h"
#include "goafem.h"
#include "mfem.hpp"
#include "mfemplus.h"
#include "parahandler.h"

#include "post.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace mfem;

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  // Check program usage
  if (argc < 2) {

    std::cout << "Usage: " << argv[0] << " input_model_filename\n";

    return 1;
  }

  ParaHandler para_handler(argv[1]);

  std::string result_path = "solutions/";
  if (access(result_path.c_str(), 0)) {
    std::cout << "folder not exists, create it ..." << std::endl;
    std::ostringstream command;
    command.str("");
    command << "mkdir " + result_path;
    int status = std::system(command.str().c_str());
  }
  Mesh *mesh = new Mesh(para_handler.mesh_file.c_str(), 1, 1);
  para_handler.preProcess();

  std::vector<DenseMatrix> cond_att(mesh->GetNE());
  for (int i = 0; i < mesh->GetNE(); i++) {
    cond_att[i] = para_handler.get_elem_conductivity(mesh->GetAttribute(i));
    mesh->SetAttribute(i, i + 1);
  }
  mesh->SetAttributes();
  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;
  GOAFEM *goafem = new GOAFEM(cond_att, para_handler, pmesh);
  // AMR loop
  for (int iter = 1; iter <= para_handler.maxit; iter++) { // iter start from 1
    double main_tic, main_toc;
    main_tic = MPI_Wtime();

    std::string amr = std::to_string(iter);
    goafem->initialize();
    if (goafem->get_problem_dofs() > para_handler.max_dofs) {
      std::cout << "\nStop due to reached max number of dofs\n";
      std::cout << "Present dofs: " << goafem->get_problem_dofs() << "\n\n";
      break;
    }

    if (goafem->get_problem_dofs() > 2147483648) {
      std::cout << "\nmax_dofs has exited the range of int-type data, please "
                   "modify your code!";
      break;
    }
    std::cout << "*********************** AMR Loop " << iter
              << "**********************\n";

    // Print number of elements, complex unknows and dofs
    goafem->print_problem_size();

    goafem->solve_primal_problem();

    std::cout << "Post processing...\n";

    Post post(pmesh, goafem->pfes, goafem->up, &para_handler);

    post.main_post_process(&para_handler, amr, goafem->sigma0,
                           goafem->cond_att);

    Meshplus pmeshplus = Meshplus(pmesh); // used to print mesh result

    goafem->solve_dual_problem();
    post.saveGridFunction(amr, "w_", &pmeshplus, goafem->w);

    goafem->error_estimating();
    Vector relative_eta = goafem->local_err;
    double max_err = relative_eta.Max();

    relative_eta *= (1.0 / max_err);

    if (para_handler.save_amr_mesh == 1) {
      post.main_save_local_mesh(amr, &pmeshplus, relative_eta, cond_att);

    } else {
      std::cout << "do not save any mesh\n" << std::endl;
    }
    goafem->refine_mesh();
    main_toc = MPI_Wtime();
    double max_time = main_toc - main_tic;
    if (iter == para_handler.maxit)
      std::cout << "the calculation time for this mesh is:\t" << max_time
                << " seconds\n";
  }
  delete goafem;
  delete pmesh;

  MPI_Finalize();

  return 0;
}
