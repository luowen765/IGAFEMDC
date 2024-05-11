

#include "goafem.h"
#include "mfemplus.h"

#include "spr.h"

#include <algorithm>
#include <cassert>
#include <cmath>

GOAFEM::GOAFEM(std::vector<DenseMatrix> &att_cond, ParaHandler &para_handler_,
               ParMesh *pmesh_)
    : cond_att(att_cond), para_handler(&para_handler_), pmesh(pmesh_),
      pfes(NULL), pfec(NULL), w(NULL), a(NULL), l(NULL), LVolume_coef(NULL),
      LGamma1_coef(NULL), integ(NULL) {
  int sn = para_handler->source_number;
  U.resize(sn);
  F.resize(sn);
  up.resize(sn);
  f.resize(sn);
  RVolume_coef.resize(sn), RGamma0_coef.resize(sn), RGamma1_coef.resize(sn);
  for (int i = 0; i < sn; i++) {
    up[i] = NULL;
    f[i] = NULL;
    RVolume_coef[i] = NULL;
    RGamma0_coef[i] = NULL;
    RGamma1_coef[i] = NULL;
  }
}

GOAFEM::~GOAFEM() {
  delete pfes;
  delete pfec;
  delete a;
  delete w;
  delete l;
  delete LVolume_coef;
  delete LGamma1_coef;
  for (int i = 0; i < para_handler->source_number; i++) {
    delete up[i];
    delete f[i];
    delete RVolume_coef[i];
    delete RGamma0_coef[i];
    delete RGamma1_coef[i];
  }
}

void GOAFEM::initialize() {

  local_sources_tets.DeleteAll();
  // find sources
  para_handler->find_point_tets(pmesh, para_handler->sources,
                                local_sources_tets, "FindPoints");
  this->set_sigma0();
  int dim = pmesh->Dimension(); // The dimention of model.
  order = 1;
  pfec = new H1_FECollection(order, dim);
  pfes = new ParFiniteElementSpace(pmesh, pfec);

  a = new ParBilinearForm(pfes);
  for (int i = 0; i < para_handler->source_number; i++) {
    // Gridfunction of primary problem
    up[i] = new ParGridFunction(pfes);
    *(up[i]) = 0.0;
    f[i] = new ParLinearForm(pfes);
    f[i]->Vector::operator=(0.0);
  }
  // Gridfunction of dual problem
  w = new ParGridFunction(pfes);
  *w = 0.0;
  // Linear form of dual problem
  l = new ParLinearForm(pfes);
  l->Vector::operator=(0.0); // VERY IMPORTANT!
  // number of elements
  local_err.SetSize(pmesh->GetNE());
  local_err = 0.0;
}

void GOAFEM::set_bdr() {
  gamma0_bdr.SetSize(pmesh->bdr_attributes.Size());
  gamma1_bdr.SetSize(pmesh->bdr_attributes.Size());
  gamma0_bdr = 0;
  gamma1_bdr = 0;
  gamma0_bdr[0] = 1;
  gamma1_bdr[1] = 1;
}

void GOAFEM::setup_integral_coefficient() {
  LVolume_coef = new LVCoefficient(cond_att);
  Vector sources_center = para_handler->get_sources_center();
  LGamma1_coef = new LGamma1Coefficient(cond_att, sources_center, pmesh);

  // Setup coefficient for right volume integral
  for (int i = 0; i < para_handler->source_number; i++) {
    Vector source_i(3);
    source_i = 0.0;
    for (int j = 0; j < 3; j++) {
      source_i[j] = para_handler->sources[i](j);
    }
    RVolume_coef[i] = new RVCoefficient(cond_att, sigma0[i], source_i);
    RGamma0_coef[i] = new RGamma0Coefficient(sigma0[i], source_i, pmesh);
    RGamma1_coef[i] =
        new RGamma1Coefficient(cond_att, sigma0[i], source_i, pmesh);
  }
}

HYPRE_Int GOAFEM::get_problem_dofs() { return pfes->GlobalTrueVSize(); }
void GOAFEM::print_problem_size() {
  HYPRE_Int size = pfes->GlobalTrueVSize();
  long n_global_elems = pmesh->GetGlobalNE();

  std::cout << "Number of tetrahedral elements: " << n_global_elems << "\n";
  std::cout << "Degrees of Freedom (DOFs): " << size << "\n";
}

void GOAFEM::solve_primal_problem() {

  std::cout << "start assemble primal problem\n";

  this->setup_integral_coefficient();
  this->set_bdr();

  integ = new DiffusionIntegrator(*LVolume_coef);
  a->AddDomainIntegrator(integ);
  a->AddBoundaryIntegrator(new MassIntegrator(*LGamma1_coef), gamma1_bdr);
  a->Assemble();
  a->Finalize();
  for (int i = 0; i < para_handler->source_number; i++) {
    f[i]->AddDomainIntegrator(new DomainLFGradIntegrator(*(RVolume_coef[i])));
    f[i]->AddBoundaryIntegrator(new BoundaryLFIntegrator(*(RGamma0_coef[i])),
                                gamma0_bdr);
    f[i]->AddBoundaryIntegrator(new BoundaryLFIntegrator(*(RGamma1_coef[i])),
                                gamma1_bdr);
    f[i]->Assemble();
    Array<int> ess_tdof_list(0);
    a->FormLinearSystem(ess_tdof_list, *(up[i]), *(f[i]), A, U[i], F[i]);
  }

  std::string problem_type = "primal";

  std::cout << "start solve primal problem\n";
  this->solve(A, U, F, problem_type);

  for (int i = 0; i < para_handler->source_number; i++) {
    a->RecoverFEMSolution(U[i], *(f[i]), *(up[i]));
  }
}

void GOAFEM::solve_dual_problem() {

  std::cout << "start assemble dual problem "
            << " \n";
  local_dual_tets.DeleteAll();

  if (para_handler->if_compute_P2 == "true") {
    std::cout << "if_compute_P2: true" << std::endl;

    for (int e = 0; e < pmesh->GetNE(); e++) {
      Element *element = pmesh->GetElement(e);
      Array<int> verts;
      element->GetVertices(verts);

      double *tt;
      for (int id = 0; id < 4; id++) {
        bool find_e = false;
        tt = pmesh->GetVertex(verts[id]);
        for (int i = 0; i < para_handler->s_plus_m.size(); i++) {
          if ((tt[0] == para_handler->s_plus_m[i](0)) &&
              (tt[1] == para_handler->s_plus_m[i](1)) &&
              (tt[2] == para_handler->s_plus_m[i](2))) {
            local_dual_tets.Append(e);
            find_e = true;
            break;
          }
        }
        if (find_e == true) {
          break;
        }
      }
    }
  } else {
    para_handler->find_point_tets(pmesh, para_handler->s_plus_m,
                                  local_dual_tets, "FindPoints");
  }

  assemble_dual_linearform();
  Array<int> ess_tdof_list(0);
  assert(a != NULL);
  // Linear system of dual problem aw=l
  a->FormLinearSystem(ess_tdof_list, *w, *l, A, W, L);

  std::string problem_type = "dual";
  std::vector<Vector> W_, L_;
  W_.push_back(W);
  L_.push_back(L);
  this->solve(A, W_, L_, problem_type);
  a->RecoverFEMSolution(W_[0], *l, *w);
}

void GOAFEM::assemble_dual_linearform() {
  // Compute non-zero element vector and assemble
  for (int i = 0; i < local_dual_tets.Size(); i++) {
    int tet_id = local_dual_tets[i];
    if (tet_id != -2) { // find it
      const FiniteElement *fe = pfes->GetFE(tet_id);
      const IntegrationRule *ir =
          &IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder() + 3);
      ElementTransformation *eltrans = pfes->GetElementTransformation(tet_id);
      double elem_volume = pmesh->GetElementVolume(tet_id);
      const int ndof = fe->GetDof();
      Vector shape(ndof);
      shape = 0.0;
      Vector elemvec(ndof);
      elemvec = 0.0;

      for (int q = 0; q < ir->GetNPoints(); ++q) {
        const IntegrationPoint &ip = ir->IntPoint(q);
        eltrans->SetIntPoint(&ip);
        double WJ = eltrans->Weight() * ip.weight;

        fe->CalcPhysShape(*eltrans, shape);
        for (int k = 0; k < ndof; k++) // loop all dofs
        {
          if (para_handler->if_use_improved_AMR == "true") {
            // add volume factor
            elemvec[k] += WJ * elem_volume * shape[k];
          } else {
            // Normalization of volume
            elemvec[k] += WJ * (1.0 / elem_volume) * shape[k];
          }
        }
      } // q loop
      Array<int> vdofs;
      pfes->GetElementDofs(tet_id, vdofs);
      l->AddElementVector(vdofs, elemvec);
    } // find it
  }
}

void GOAFEM::solve_with_pcg(OperatorHandle &A, Vector &X, Vector &B,
                            std::string pre, std::string problem_type) {
  double tol = 0.0;
  if (problem_type == "primal") {
    tol = 1e-8;
  } else if (problem_type == "dual") {
    tol = 1e-4;
  } else {
    std::cout << "please check your problem_type! primal or dual?" << std::endl;
    std::abort();
  }

  Solver *prec = NULL;

  if (pre == "ilu") {
    HypreILU *ilu = new HypreILU;
    ilu->SetPrintLevel(0);
    prec = ilu;
  } else {
    std::cout << "Unsupported input linear solver: "
              << para_handler->linear_solver << "\n";
    std::abort();
  }

  CGSolver cg(MPI_COMM_WORLD);
  cg.SetMaxIter(1000);
  cg.SetRelTol(tol);
  cg.SetPrintLevel(0);

  if (prec) {
    cg.SetPreconditioner(*prec);
  }
  cg.SetOperator(*A);
  cg.Mult(B, X);

  delete prec;
}

void GOAFEM::solve(OperatorHandle &A, std::vector<Vector> &X,
                   std::vector<Vector> &B, std::string problem_type) {

  std::string linear_solver = para_handler->linear_solver;
  if (linear_solver == "mumps") {
    MUMPSSolver mumps(MPI_COMM_WORLD);
    mumps.SetPrintLevel(0);
    mumps.SetMatrixSymType(MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
    mumps.SetOperator(*A.Ptr());
    for (int i = 0; i < B.size(); i++) {
      mumps.Mult(B[i], X[i]);
    }
  } else {
    for (int i = 0; i < B.size(); i++) {
      solve_with_pcg(A, X[i], B[i], linear_solver, problem_type);
    }
  }
}

void GOAFEM::set_sigma0() {

  int ns = local_sources_tets.Size();
  assert(ns == para_handler->source_number);
  global_initial_sources_tets.SetSize(ns);

  // if (myid == 0) {
  for (int i = 0; i < ns; i++) {
    int tet_id = local_sources_tets[i];
    if (tet_id != -2) {
      global_initial_sources_tets[i] = pmesh->GetAttribute(tet_id) - 1;
    } else {
      std::cout << "some source points are not found!\n";
      std::abort();
    }
  }

  sigma0.clear();
  sigma0.resize(para_handler->source_number);

  for (int k = 0; k < para_handler->source_number; k++) {
    sigma0[k] = cond_att[global_initial_sources_tets[k]];
  }
}

void GOAFEM::error_estimating() {

  assert(local_err.Size() == pmesh->GetNE());
  int c = para_handler->source_number;
  int dofs = pmesh->GetNV();

  SPR spr(para_handler, pmesh, pfes);
  std::vector<DenseVector> Ui(c);

  for (int i = 0; i < c; i++) {
    Ui[i].resize(dofs);
    for (int j = 0; j < dofs; j++) {
      Ui[i](j) = (*up[i])(j);
    }
  }
  DenseVector W(dofs);
  for (int j = 0; j < dofs; j++) {

    W[j] = (*w)(j);
  }

  spr.gradient_patch_recovery(Ui, W, cond_att, local_err);
  std::cout << "errors.Max(): " << local_err.Max() << std::endl;
  std::cout << "errors.mean(): " << local_err.Sum() / local_err.Size()
            << std::endl;
}

void GOAFEM::refine_mesh() {

  double threshold;

  if (para_handler->amr_type == 0) {
    int N_cell = pmesh->GetNE();
    std::vector<double> eta(N_cell);
    for (int i = 0; i < N_cell; i++) {
      eta[i] = local_err[i];
    }
    sort(std::begin(eta), std::end(eta));
    threshold = eta[round((1 - para_handler->beta) * N_cell)];
  } else {
    threshold = para_handler->beta * local_err.Max();
  }

  pmesh->RefineByError(local_err, threshold);
}
