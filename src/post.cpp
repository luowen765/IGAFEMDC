
#include "post.h"
#include "dc.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
// #include <mpi.h>
#include <unistd.h>

Post::Post(ParMesh *pmesh_, ParFiniteElementSpace *pfes_,
           std::vector<ParGridFunction *> &up_, ParaHandler *para_handler_)
    : pmesh(pmesh_), pfes(pfes_), up(up_), para_handler(para_handler_) {}

Post::~Post() {}

// Post-processing and save the solution
void Post::main_post_process(ParaHandler *para_handler, std::string amr,
                             std::vector<DenseMatrix> &sigma0,
                             std::vector<DenseMatrix> &cond_att) {

  local_sites_tets.resize(para_handler->source_number);

  if (para_handler->if_compute_P2 == "true") {

    this->compute_oneSource_U_E_J_averagely(cond_att, sigma0);

  } else {

    local_sites_tets.resize(para_handler->source_number);

    for (int i = 0; i < para_handler->source_number; i++) {

      // Find elements which contain sites
      para_handler->find_point_tets(pmesh, para_handler->sites[i],
                                    local_sites_tets[i], "FindPoints");
    }

    // compute Up of measurements
    this->post_process(cond_att, sigma0);
  }

  std::string output_file = "solutions/solution" + std::string(".") + amr;
  // this->read_compute_survery_mode(para_handler->sorted_survey, output_file);
  this->read_compute_survery_mode(para_handler->survey_input, output_file);
}

void Post::post_process(std::vector<DenseMatrix> &cond_att,
                        std::vector<DenseMatrix> &sigma0) {

  int Ns = para_handler->source_number;
  local_Up.resize(Ns);

  if (para_handler->if_compute_J == "true") {
    Ex.resize(Ns);
    Ey.resize(Ns);
    Ez.resize(Ns);
    Jx.resize(Ns);
    Jy.resize(Ns);
    Jz.resize(Ns);
  }

  for (int k = 0; k < para_handler->source_number; k++) {
    int nsite = para_handler->sites[k].size();
    assert(sigma0.size() == para_handler->source_number);
    local_Up[k].SetSize(nsite);
    local_Up[k] = 0.0;

    if (para_handler->if_compute_J == "true") {
      Ex[k].SetSize(nsite);
      Ey[k].SetSize(nsite);
      Ez[k].SetSize(nsite);
      Jx[k].SetSize(nsite);
      Jy[k].SetSize(nsite);
      Jz[k].SetSize(nsite);
    }

    for (int i = 0; i < nsite; i++) {
      int tet_id = local_sites_tets[k][i];
      if (tet_id != -2) { // find in this process

        const FiniteElement *fe = pfes->GetFE(tet_id);
        int ndof = fe->GetDof();
        Vector pt(3);
        pt = 0.0;
        pt[0] = para_handler->sites[k][i](0);
        pt[1] = para_handler->sites[k][i](1);
        pt[2] = para_handler->sites[k][i](2);
        IntegrationPoint ip;
        ElementTransformation *tr = pfes->GetElementTransformation(tet_id);
        tr->TransformBack(pt, ip);
        tr->SetIntPoint(&ip);
        // Compute shape function
        Vector shape(ndof);
        shape = 0.0;
        fe->CalcPhysShape(*tr, shape);
        Array<int> vdofs;
        pfes->GetElementDofs(tet_id, vdofs);
        Vector tempUp;
        tempUp = 0.0;
        up[k]->GetSubVector(vdofs, tempUp);

        double potential = 0.0;
        for (int j = 0; j < ndof; j++) {
          potential += shape[j] * tempUp[j];
        }
        local_Up[k][i] = potential; // local Up

        if (para_handler->if_compute_J == "true") {
          DenseMatrix Dshape(ndof, 3); // 4*3
          Dshape = 0.0;
          fe->CalcPhysDShape(*tr, Dshape);
          Vector tempE(3);
          tempE = 0.0;
          Dshape.AddMultTranspose(tempUp, tempE);
          Vector deltaUs(3);
          deltaUs = 0.0;
          Vector pi(3);
          pi = 0.0;
          for (int kk = 0; kk < 3; kk++) {
            pi[kk] = para_handler->sources[k](kk);
          }
          if (equalVector(pt, pi)) {
            std::cout << "pt: " << std::endl;
            pt.Print();
            std::cout << "pi: " << std::endl;
            pi.Print();
          }
          assert(!equalVector(pt, pi));
          gradient_U_i_s(pt, pi, sigma0[k], deltaUs);
          assert(Ex[k].Size() != 0);
          Ex[k][i] = -1.0 * (deltaUs(0) + tempE(0));
          Ey[k][i] = -1.0 * (deltaUs(1) + tempE(1));
          Ez[k][i] = -1.0 * (deltaUs(2) + tempE(2));
          Vector totalE(3), totalJ(3);
          totalE = 0.0;
          totalJ = 0.0;
          totalE[0] = Ex[k][i];
          totalE[1] = Ey[k][i];
          totalE[2] = Ez[k][i];
          DenseMatrix tet_cond = cond_att[tr->Attribute - 1];
          tet_cond.AddMult(totalE, totalJ);
          Jx[k][i] = totalJ(0);
          Jy[k][i] = totalJ(1);
          Jz[k][i] = totalJ(2);
        }

      } else {
        std::cout << "some points are not found!\n";
        std::abort();
      }
    }
  }

  Ut.resize(Ns);
  for (int k = 0; k < Ns; k++) {
    int sites_number = para_handler->sites[k].size();
    Ut[k].SetSize(sites_number);
    Vector pi(3);
    pi = 0.0;
    for (int kk = 0; kk < 3; kk++) {
      pi[kk] = para_handler->sources[k](kk);
    }
    for (int h = 0; h < sites_number; h++) {
      Vector pj(3);
      pj = 0.0;
      for (int gg = 0; gg < 3; gg++) {
        pj[gg] = para_handler->sites[k][h](gg);
      }
      if (pj(0) == pi(0) && pj(1) == pi(1) && pj(2) == pi(2)) {
        Ut[k][h] = 1E10;
      } else {
        Ut[k][h] = local_Up[k][h] + U_i_s(pj, pi, sigma0[k]);
      }
    }

  } // source loop
}

// Save local mesh
void Post::main_save_local_mesh(std::string amr_, Meshplus *pmesh,
                                Vector &local_eta,
                                std::vector<DenseMatrix> &att_cond) {

  MPI_Barrier(MPI_COMM_WORLD);
  std::ostringstream vtk_name2;
  vtk_name2 << "solutions/"
            << "err" << std::setfill('0') << std::setw(2) << "." << amr_
            << ".vtk";
  std::ofstream vtk_ofs2(vtk_name2.str().c_str());
  vtk_ofs2.precision(8);
  pmesh->PrintVTK_eta(vtk_ofs2, local_eta);
}
// checked
void Post::read_compute_survery_mode(std::string survery_mode_file,
                                     std::string output_file) {

  unsigned int n_data = 0;
  unsigned int mode = 0;
  int source;

  std::ifstream mode_stream(survery_mode_file.c_str());
  assert(mode_stream.good());

  std::ofstream out_stream(output_file.c_str());
  assert(out_stream.good());

  mode_stream >> n_data >> mode;

  if (mode == 11) { // pole-pole configuration
    out_stream << n_data << "\t" << mode << "\n";
    int ns = para_handler->sources.size();
    int surveys = 0;
    for (int i = 0; i < ns; i++) {
      int nsite = para_handler->sites[i].size();
      surveys += nsite;
      for (int j = 0; j < nsite; j++) {
        int C_local_id = i;
        int P_local_id = j;

        Vertex C_node = para_handler->sources[C_local_id];
        Vertex P_node = para_handler->sites[C_local_id][P_local_id];
        double a = length_two_point(C_node.operator()(), P_node.operator()());
        double k = 2.0 * PI * a;
        double U = this->Ut[C_local_id][P_local_id];
        double rho = k * U;
        out_stream << C_local_id << "\t" << (C_node)(0) << "\t" << (C_node)(1)
                   << "\t" << (C_node)(2) << "\t" << P_local_id << "\t"
                   << (P_node)(0) << "\t" << (P_node)(1) << "\t" << (P_node)(2)
                   << "\t" << k << "\t" << U << "\t" << rho << "\n";
      }
    }
    assert(n_data == surveys);

  }

  else if (mode = 41) {

    if (para_handler->if_compute_P2 == "true") {

      out_stream << "s1x\t"
                 << "s1y\t"
                 << "s1z\t"
                 << "s2x\t"
                 << "s2y\t"
                 << "s2z\t"
                 << "s3x\t"
                 << "s3y\t"
                 << "s3z\t"
                 << "s4x\t"
                 << "s4y\t"
                 << "s4z\t"
                 << "x\t"
                 << "y\t"
                 << "z\t"
                 << "U\t"
                 << "rho1\t"
                 << "rho2\t"
                 << "P2\n";
      int ns = para_handler->sources.size() / 4;
      int surveys = 0;
      for (int i = 0; i < ns; i++) {
        int nsite = para_handler->sites[i].size();
        surveys += nsite;
        for (int j = 0; j < nsite; j++) {
          int c1_local = i;
          int c2_local = i + 1;
          int c3_local = i + 2;
          int c4_local = i + 3;
          int p_local = j;
          Vertex C1_node = para_handler->sources[c1_local];
          Vertex C2_node = para_handler->sources[c2_local];
          Vertex C3_node = para_handler->sources[c3_local];
          Vertex C4_node = para_handler->sources[c4_local];
          Vertex P_node = para_handler->sites[c1_local][p_local];
          // compute Ut ----------------
          double Ut_s1 = Ut[c1_local][p_local] - Ut[c2_local][p_local];
          double Ut_s2 = Ut[c3_local][p_local] - Ut[c4_local][p_local];
          double Ut_4source = Ut_s1 + Ut_s2;

          //-----------------------
          // compute rho
          double AM =
              length_two_point(C1_node.operator()(), P_node.operator()());
          double BM =
              length_two_point(C2_node.operator()(), P_node.operator()());
          double CM =
              length_two_point(C3_node.operator()(), P_node.operator()());
          double DM =
              length_two_point(C4_node.operator()(), P_node.operator()());
          double k1 = 2.0 * PI / (1.0 / AM - 1.0 / BM);
          double k2 = 2.0 * PI / (1.0 / CM - 1.0 / DM);
          double rho1 = k1 * Ut_s1;
          double rho2 = k2 * Ut_s2;

          double P2 = -1.0;
          // compute P2 variable ----------------
          assert(Ex[0].Size() != 0);
          double E11 = Ex[c1_local][p_local] - Ex[c2_local][p_local];
          double E12 = Ey[c1_local][p_local] - Ey[c2_local][p_local];
          double E21 = Ex[c3_local][p_local] - Ex[c4_local][p_local];
          double E22 = Ey[c3_local][p_local] - Ey[c4_local][p_local];
          double RA = sqrt(pow(C1_node(0) - P_node(0), 2) +
                           pow(C1_node(1) - P_node(1), 2));
          double RB = sqrt(pow(C2_node(0) - P_node(0), 2) +
                           pow(C2_node(1) - P_node(1), 2));
          double RC = sqrt(pow(C3_node(0) - P_node(0), 2) +
                           pow(C3_node(1) - P_node(1), 2));
          double RD = sqrt(pow(C4_node(0) - P_node(0), 2) +
                           pow(C4_node(1) - P_node(1), 2));
          double J11 = ((P_node(0) - C1_node(0)) * (1.0 / (RA * RA * RA)) -
                        (P_node(0) - C2_node(0)) * (1.0 / (RB * RB * RB))) *
                       (1.0 / (2.0 * PI));
          double J12 = ((P_node(1) - C1_node(1)) * (1.0 / (RA * RA * RA)) -
                        (P_node(1) - C2_node(1)) * (1.0 / (RB * RB * RB))) *
                       (1.0 / (2.0 * PI));
          double J21 = ((P_node(0) - C3_node(0)) * (1.0 / (RC * RC * RC)) -
                        (P_node(0) - C4_node(0)) * (1.0 / (RD * RD * RD))) *
                       (1.0 / (2.0 * PI));
          double J22 = ((P_node(1) - C3_node(1)) * (1.0 / (RC * RC * RC)) -
                        (P_node(1) - C4_node(1)) * (1.0 / (RD * RD * RD))) *
                       (1.0 / (2.0 * PI));
          double J_det = std::abs(J11 * J22 - J21 * J12);
          double rho11 = (E11 * J22 - E21 * J12) / J_det;
          double rho12 = (E21 * J11 - E11 * J21) / J_det;
          double rho21 = (E12 * J22 - E22 * J12) / J_det;
          double rho22 = (E22 * J11 - E12 * J21) / J_det;

          P2 = sqrt(abs(rho11 * rho22 - rho12 * rho21));

          //----------------------------------------
          out_stream << (C1_node)(0) << "\t" << (C1_node)(1) << "\t"
                     << (C1_node)(2) << "\t" << (C2_node)(0) << "\t"
                     << (C2_node)(1) << "\t" << (C2_node)(2) << "\t"
                     << (C3_node)(0) << "\t" << (C3_node)(1) << "\t"
                     << (C3_node)(2) << "\t" << (C4_node)(0) << "\t"
                     << (C4_node)(1) << "\t" << (C4_node)(2) << "\t"
                     << (P_node)(0) << "\t" << (P_node)(1) << "\t"
                     << (P_node)(2) << "\t" << Ut_4source << "\t" << rho1
                     << "\t" << rho2 << "\t" << P2 << "\n";
        }
      }
      assert(n_data == surveys);
    }

    if (para_handler->if_compute_J == "true") {

      std::string out_JU = "EJU.txt";

      std::ofstream out(out_JU);
      assert(out.good());
      out.precision(8);
      out << "x"
          << "\t"
          << "y"
          << "\t"
          << "z"
          << "\t"
          << "Ex"
          << "\t"
          << "Ey"
          << "\t"
          << "Ez"
          << "\t"
          << "Jx"
          << "\t"
          << "Jy"
          << "\t"
          << "Jz"
          << "\t"
          << "U\n";

      int Nsite = para_handler->sites[0].size();
      for (int j = 0; j < Nsite; j++) {
        out << para_handler->sites[0][j](0) << "\t"
            << para_handler->sites[0][j](1) << "\t"
            << para_handler->sites[0][j](2) << "\t"
            << Ex[0][j] - Ex[1][j] + Ex[2][j] - Ex[3][j] << "\t"
            << Ey[0][j] - Ey[1][j] + Ey[2][j] - Ey[3][j] << "\t"
            << Ez[0][j] - Ez[1][j] + Ez[2][j] - Ez[3][j] << "\t"
            << Jx[0][j] - Jx[1][j] + Jx[2][j] - Jx[3][j] << "\t"
            << Jy[0][j] - Jy[1][j] + Jy[2][j] - Jy[3][j] << "\t"
            << Jz[0][j] - Jz[1][j] + Jz[2][j] - Jz[3][j] << "\t"
            << Ut[0][j] - Ut[1][j] + Ut[2][j] - Ut[3][j]
            //  <<"\t"<<local_Up[i][j]
            << std::endl;
      }
    }
  }

  else {
    std::cout << "Other modes are not supported temporarily!";
  }

  mode_stream.close();
  out_stream.close();

  return;
}

void Post::compute_oneSource_U_E_J_averagely(std::vector<DenseMatrix> &cond_att,
                                             std::vector<DenseMatrix> &sigma0) {
  assert(sigma0.size() == para_handler->source_number);
  int Ns = para_handler->source_number;

  std::vector<std::map<Vertex *, std::vector<double>>> data(Ns);
  std::vector<std::map<Vertex *, std::vector<double>>> data_ex(Ns);
  std::vector<std::map<Vertex *, std::vector<double>>> data_ey(Ns);
  std::vector<std::map<Vertex *, std::vector<double>>> data_ez(Ns);
  std::vector<std::map<Vertex *, std::vector<double>>> data_jx(Ns);
  std::vector<std::map<Vertex *, std::vector<double>>> data_jy(Ns);
  std::vector<std::map<Vertex *, std::vector<double>>> data_jz(Ns);

  for (int e = 0; e < pmesh->GetNE(); e++) {
    Element *element = pmesh->GetElement(e);
    Array<int> verts;
    element->GetVertices(verts);

    double *tt;
    for (int id = 0; id < 4; id++) {
      bool find_e = false;
      tt = pmesh->GetVertex(verts[id]);
      for (int ic = 0; ic < Ns; ic++) {
        for (int im = 0; im < para_handler->sites[ic].size(); im++) {
          if ((tt[0] == para_handler->sites[ic][im](0)) &&
              (tt[1] == para_handler->sites[ic][im](1)) &&
              (tt[2] == para_handler->sites[ic][im](2))) {
            Vertex *node = NULL;
            node = &(para_handler->sites[ic][im]);
            const FiniteElement *fe = pfes->GetFE(e);
            int ndof = fe->GetDof();

            Vector pt(3);
            pt = 0.0;
            pt[0] = (*node)(0);
            pt[1] = (*node)(1);
            pt[2] = (*node)(2);
            Vector pi(3);
            pi = 0.0;
            for (int kk = 0; kk < 3; kk++) {
              pi[kk] = para_handler->sources[ic](kk);
            }
            IntegrationPoint ip;
            ElementTransformation *tr = pfes->GetElementTransformation(e);
            tr->TransformBack(pt, ip);
            tr->SetIntPoint(&ip);

            Array<int> vdofs;

            pfes->GetElementDofs(e, vdofs);
            Vector tempUp;
            tempUp = 0.0;
            up[ic]->GetSubVector(vdofs, tempUp);

            Vector shape(ndof);
            shape = 0.0;
            fe->CalcPhysShape(*tr, shape);

            double potential = 0.0;
            for (int j = 0; j < ndof; j++) {
              potential += shape[j] * tempUp[j];
            }

            double total_potential = 0.0;
            if (pt(0) == pi(0) && pt(1) == pi(1) && pt(2) == pi(2)) {
              total_potential = 1E10;
            } else {
              total_potential = potential + U_i_s(pt, pi, sigma0[ic]);
            }

            // compute E and J
            Vector tempE(3);
            tempE = 0.0;
            DenseMatrix Dshape(ndof, 3); // 4*3
            Dshape = 0.0;
            fe->CalcPhysDShape(*tr, Dshape);
            Dshape.AddMultTranspose(tempUp, tempE);
            // tempE *= -1.0;
            DenseMatrix tet_cond = cond_att[tr->Attribute - 1];
            Vector deltaUs(3);
            deltaUs = 0.0;
            if (equalVector(pt, pi)) {
              std::cout << "pt: " << std::endl;
              pt.Print();
              std::cout << "pi: " << std::endl;
              pi.Print();
            }
            assert(!equalVector(pt, pi));

            gradient_U_i_s(pt, pi, sigma0[ic], deltaUs);
            Vector totalE(3), totalJ(3);
            totalE = 0.0;
            totalJ = 0.0;
            totalE[0] = (tempE(0) + deltaUs(0)) * -1.0;
            totalE[1] = (tempE(1) + deltaUs(1)) * -1.0;
            totalE[2] = (tempE(2) + deltaUs(2)) * -1.0;
            tet_cond.AddMult(totalE, totalJ);

            data[ic][node].push_back(total_potential);
            data_ex[ic][node].push_back(totalE(0));
            data_ey[ic][node].push_back(totalE(1));
            data_ez[ic][node].push_back(totalE(2));
            data_jx[ic][node].push_back(totalJ(0));
            data_jy[ic][node].push_back(totalJ(1));
            data_jz[ic][node].push_back(totalJ(2));

            break;
          }
        }
      }
    }
  }

  Ut.resize(Ns);

  Ex.resize(Ns);
  Ey.resize(Ns);
  Ez.resize(Ns);
  Jx.resize(Ns);
  Jy.resize(Ns);
  Jz.resize(Ns);

  for (int ic = 0; ic < para_handler->source_number; ic++) {
    int nsite = para_handler->sites[ic].size();
    for (int im = 0; im < nsite; im++) {
      Ut[ic].SetSize(nsite);

      Ex[ic].SetSize(nsite);
      Ey[ic].SetSize(nsite);
      Ez[ic].SetSize(nsite);
      Jx[ic].SetSize(nsite);
      Jy[ic].SetSize(nsite);
      Jz[ic].SetSize(nsite);
    }
    typedef std::map<Vertex *, std::vector<double>>::iterator IT;
    unsigned int p = 0;
    unsigned int pex = 0;
    unsigned int pey = 0;
    unsigned int pez = 0;
    unsigned int pjx = 0;
    unsigned int pjy = 0;
    unsigned int pjz = 0;
    std::map<Vertex *, std::vector<double>> &node2potential = data[ic];
    // map内元素自动排序
    for (IT it = node2potential.begin(); it != node2potential.end();
         it++, p++) {
      std::vector<double> &temp = (*it).second;
      double summation = std::accumulate(temp.begin(), temp.end(), 0.0);
      double averaged_potential = summation / temp.size();
      Ut[ic][p] = averaged_potential;
    }

    std::map<Vertex *, std::vector<double>> &node2potential_ex = data_ex[ic];
    std::map<Vertex *, std::vector<double>> &node2potential_ey = data_ey[ic];
    std::map<Vertex *, std::vector<double>> &node2potential_ez = data_ez[ic];
    std::map<Vertex *, std::vector<double>> &node2potential_jx = data_jx[ic];
    std::map<Vertex *, std::vector<double>> &node2potential_jy = data_jy[ic];
    std::map<Vertex *, std::vector<double>> &node2potential_jz = data_jz[ic];
    for (IT it_ex = node2potential_ex.begin(); it_ex != node2potential_ex.end();
         it_ex++, pex++) {
      std::vector<double> &temp_ex = (*it_ex).second;
      double summation_ex =
          std::accumulate(temp_ex.begin(), temp_ex.end(), 0.0);
      double averaged_potential_ex = summation_ex / temp_ex.size();
      Ex[ic][pex] = averaged_potential_ex;
    }

    for (IT it_ey = node2potential_ey.begin(); it_ey != node2potential_ey.end();
         it_ey++, pey++) {
      std::vector<double> &temp_ey = (*it_ey).second;
      double summation_ey =
          std::accumulate(temp_ey.begin(), temp_ey.end(), 0.0);
      double averaged_potential_ey = summation_ey / temp_ey.size();
      Ey[ic][pey] = averaged_potential_ey;
    }

    for (IT it_ez = node2potential_ez.begin(); it_ez != node2potential_ez.end();
         it_ez++, pez++) {
      std::vector<double> &temp_ez = (*it_ez).second;
      double summation_ez =
          std::accumulate(temp_ez.begin(), temp_ez.end(), 0.0);
      double averaged_potential_ez = summation_ez / temp_ez.size();
      Ez[ic][pez] = averaged_potential_ez;
    }

    for (IT it_jx = node2potential_jx.begin(); it_jx != node2potential_jx.end();
         it_jx++, pjx++) {
      std::vector<double> &temp_jx = (*it_jx).second;
      double summation_jx =
          std::accumulate(temp_jx.begin(), temp_jx.end(), 0.0);
      double averaged_potential_jx = summation_jx / temp_jx.size();
      Jx[ic][pjx] = averaged_potential_jx;
    }

    for (IT it_jy = node2potential_jy.begin(); it_jy != node2potential_jy.end();
         it_jy++, pjy++) {
      std::vector<double> &temp_jy = (*it_jy).second;
      double summation_jy =
          std::accumulate(temp_jy.begin(), temp_jy.end(), 0.0);
      double averaged_potential_jy = summation_jy / temp_jy.size();
      Jy[ic][pjy] = averaged_potential_jy;
    }

    for (IT it_jz = node2potential_jz.begin(); it_jz != node2potential_jz.end();
         it_jz++, pjz++) {
      std::vector<double> &temp_jz = (*it_jz).second;
      double summation_jz =
          std::accumulate(temp_jz.begin(), temp_jz.end(), 0.0);
      double averaged_potential_jz = summation_jz / temp_jz.size();
      Jz[ic][pjz] = averaged_potential_jz;
    }
  }
}

void Post::saveGridFunction(std::string amr_, std::string fn, Meshplus *pmesh,
                            ParGridFunction *field) {

  std::ostringstream vtk_name;
  vtk_name << "solutions/" << fn << amr_ << ".vtk";

  std::ofstream vtk_ofs(vtk_name.str().c_str());
  vtk_ofs.precision(8);

  // field->SaveVTK(vtk_ofs, "w",1);
  pmesh->write_out_node_info_vtk(vtk_ofs, *field);
}