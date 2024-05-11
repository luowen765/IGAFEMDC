// C++ include
#include <time.h>
// Local include
#include "spr.h"
#include <iostream>

//--------------------------------------------------------------------------
SPR::SPR(ParaHandler *para_handler_, ParMesh *pmesh,
         ParFiniteElementSpace *pfes, const unsigned int s)
    : _para_handler(para_handler_), _pmesh(pmesh), _pfes(pfes), _s(s) {}

void SPR::gradient_patch_recovery(std::vector<DenseVector> &Ui, DenseVector &W,
                                  std::vector<DenseMatrix> &cond_att,
                                  Vector &local_err) {

  int Ne = _pmesh->GetNE();

  Ui.push_back(W); // c+1 solution
  std::set<unsigned int> computed_elems;

  for (unsigned int i = 0; i < Ne; i++) {

    unsigned int id = i;

    if (computed_elems.find(id) == computed_elems.end()) { // already estimated
      // Build the element patch surrounding "elem"
      std::set<unsigned int> patch_id;

      this->get_element_patch(id, patch_id, this->_s, cond_att);

      assert(patch.size() > 0); // for linear.
      // C is a 4x4 symmetrical matrix for linear polynomials over path
      Matrix4D C;
      std::vector<Matrix43D> b, a;

      this->compute_C(patch_id, C);
      // std::cout << "C: " << C << std::endl;
      this->get_patch_b(patch_id, Ui, b);
      // std::cout << "b: \n" << b[0] << std::endl;
      this->solve_a(C, b, a);

      std::map<unsigned int, double> energy_norm;
      this->compute_PI_k(patch_id, Ui, a, cond_att, energy_norm);

      std::set<unsigned int>::const_iterator it = patch_id.begin();
      for (; it != patch_id.end(); it++) {

        local_err[(*it)] = energy_norm[(*it)];
        computed_elems.insert((*it));
      }
    }
  }

  Ui.erase(Ui.end());
  std::vector<double> temp(W.size());
  for (int i = 0; i < temp.size(); i++)
    temp[i] = W(i);

  return;
}

// checked
void SPR::get_element_patch(unsigned int elem,
                            // const unsigned int elem,
                            std::set<unsigned int> &patch_id,
                            const unsigned int S,
                            std::vector<DenseMatrix> &cond_att) {
  assert(elem != NULL);
  // First add the element of itself
  patch_id.clear();
  patch_id.insert(elem);
  const unsigned int N = 20;
  // const unsigned int N = 2;
  std::set<unsigned int> neighbors;
  neighbors.insert(elem);

  for (int i = 0; i < N; i++) {
    // Loop over neighbors
    std::set<unsigned int>::const_iterator it = neighbors.begin();
    const std::set<unsigned int>::const_iterator end = neighbors.end();
    std::set<unsigned int> temp;
    for (; it != end; ++it) {
      int elem_temp_id = *it;
      Array<int> faces, ori;
      _pmesh->GetElementFaces(elem_temp_id, faces, ori);
      for (unsigned int s = 0; s < faces.Size(); s++) {
        int face_id = faces[s];
        int nr_id = -1;
        bool haveNeighbor = false;
        if (_pmesh->FaceIsInterior(face_id)) {
          haveNeighbor = true;
          int elem1_id, elem2_id;
          _pmesh->GetFaceElements(face_id, &elem1_id, &elem2_id);
          assert(elem1_id != -1 && elem2_id != -1);
          if (elem_temp_id == elem1_id) {
            nr_id = elem2_id;
          } else {
            assert(elem_temp_id == elem2_id);
            nr_id = elem1_id;
          }
        }
        if (haveNeighbor) { // we have a neighbor on this side
          const Element *neighbor = _pmesh->GetElement(nr_id);
          const Element *e = _pmesh->GetElement(elem);

          DenseMatrix temp1(3, 3), temp2(3, 3);
          temp1 = cond_att[neighbor->GetAttribute() - 1];
          temp2 = cond_att[e->GetAttribute() - 1];
          if (equalDM(temp1, temp2)) { 
            patch_id.insert(nr_id);
            temp.insert(nr_id);
          }
        }
      }
    }
    neighbors = temp;
    if (patch_id.size() > S)
      break;
  } // end while loop
  // Waring if patch.size()<4
  if (patch_id.size() < 4)
    std::cout << "Waring,in SPR::get_element_patch()"
              << "\nThe size of patch is < 4\n";
}

void SPR::get_patch_polynomial_c(std::vector<double> &c, const Point_ren &p) {
  const double x = p(0), y = p(1), z = p(2);
  c.clear();
  c.push_back(1); // a0
  c.push_back(x); // a1
  c.push_back(y); // a2
  c.push_back(z); // a2
  return;
}

void SPR::compute_C(std::set<unsigned int> &patch_id, Matrix4D &C) {
  assert(patch.size() >= 4);
  C.setZero();

  // Loop the elements in patch.
  std::set<unsigned int>::const_iterator it = patch_id.begin();
  for (; it != patch_id.end(); it++) {
    int elem_id = (*it);

    // int elem_id = (*it);
    ElementTransformation *eltrans = _pfes->GetElementTransformation(elem_id);
    const FiniteElement *fe = _pfes->GetFE(elem_id);
    const IntegrationRule *ir =
        &IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder());
    for (unsigned int q = 0; q < ir->GetNPoints(); q++) {

      const IntegrationPoint &ip = ir->IntPoint(q);
      eltrans->SetIntPoint(&ip); // eltrans不能是const值
      double WJ = eltrans->Weight() * ip.weight;

      std::vector<double> c;
      Vector r_temp(3);
      r_temp = 0.0;
      eltrans->Transform(ip, r_temp);
      Point_ren r(r_temp[0], r_temp[1], r_temp[2]);
      // Get the value of "c" at Gauss Point_ren "r"
      this->get_patch_polynomial_c(c, r);
      // Check the size
      assert(c.size() == 4);
      // Calculate the "C" in equation (12) of our GEOPHYSICS paper.
      for (unsigned int i = 0; i < 4; i++)
        for (unsigned int j = 0; j < 4; j++)
          C(i, j) += WJ * c[i] * c[j];
    }
  }

  return;
}

void SPR::get_patch_b(std::set<unsigned int> &patch_id,
                      std::vector<DenseVector> &solution,
                      std::vector<Matrix43D> &b) {
  // The a is got by solving the small equation Ca=b, where C is 4*4 matrix
  // and b are 4*3 matrix.
  assert(patch.size() >= 4);
  unsigned int n_solution = solution.size();
  b.resize(n_solution);

  // Loop the elements in patch.
  std::set<unsigned int>::const_iterator it = patch_id.begin();
  // std::cout << "patch.size(): " << patch.size() << std::endl;
  for (; it != patch_id.end(); it++) {
    int elem_id = (*it);
    ElementTransformation *eltrans = _pfes->GetElementTransformation(elem_id);
    const FiniteElement *fe = _pfes->GetFE(elem_id);
    const IntegrationRule *ir =
        &IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder());

    std::vector<unsigned int> DOFs;

    int ndof = fe->GetDof();
    int dim = fe->GetDim();
    Array<int> vdofs;
    _pfes->GetElementDofs(elem_id, vdofs);

    DOFs.resize(vdofs.Size());
    for (int i = 0; i < vdofs.Size(); i++) {
      DOFs[i] = vdofs[i];
    }
    assert(DOFs.size() == 4);

    for (unsigned int q = 0; q < ir->GetNPoints(); q++) {

      const IntegrationPoint &ip = ir->IntPoint(q);
      eltrans->SetIntPoint(&ip);

      double WJ = eltrans->Weight() * ip.weight;

      DenseMatrix Dshape(ndof, 3);
      Dshape = 0.0;
      fe->CalcPhysDShape(*eltrans, Dshape);

      Point_ren dev_shape[4];
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
          dev_shape[i](j) = Dshape(i, j);
        }
      }

      Vector r_temp(3);
      r_temp = 0.0;
      eltrans->Transform(ip, r_temp);
      Point_ren r(r_temp[0], r_temp[1], r_temp[2]);

      std::vector<double> c;
      // Get the value of "c" at Gauss poin-t "r"
      this->get_patch_polynomial_c(c, r);

      for (int n = 0; n < n_solution; n++) {
        Point_ren grad_u_h;
        b[n].setZero();
        // Calculate the gradients at gauss poin-ts.
        for (unsigned int i = 0; i < 4; i++)
          grad_u_h = grad_u_h + dev_shape[i] * solution[n](DOFs[i]);

        for (unsigned int i = 0; i < 4; i++)
          for (unsigned int j = 0; j < 3; j++) {
            b[n](i, j) += WJ * c[i] * grad_u_h(j);
          }
      }
    }
  }

  return;
}

void SPR::solve_a(Matrix4D &C, std::vector<Matrix43D> &b,
                  std::vector<Matrix43D> &a) {
  a.resize(b.size());
  for (int i = 0; i < a.size(); i++)
    a[i] = C.partialPivLu().solve(b[i]);

  return;
}

void SPR::compute_P_value(Point_ren &r, Matrix43D &a, Matrix13D &P) {
  P.setZero();
  std::vector<double> c;
  this->get_patch_polynomial_c(c, r);
  assert(c.size() == 4);
  Matrix14D dv_c;
  for (int t = 0; t < 4; t++)
    dv_c(t) = c[t];
  P = dv_c * a;

  return;
}

void SPR::compute_PI_k(std::set<unsigned int> &patch_id,
                       std::vector<DenseVector> &solution,
                       std::vector<Matrix43D> &a,
                       //  std::vector<Matrix3D> &sigma,
                       std::vector<DenseMatrix> &cond_att,
                       std::map<unsigned int, double> &energy_norm) {
  assert(solution.size() == a.size());
  energy_norm.clear();
  // solution is Ui, W

  // Loop the elements in patch.
  std::set<unsigned int>::const_iterator it = patch_id.begin();
  for (; it != patch_id.end(); it++) {
    int elem_id = (*it);
    ElementTransformation *eltrans = _pfes->GetElementTransformation(elem_id);
    const FiniteElement *fe = _pfes->GetFE(elem_id);
    const IntegrationRule *ir =
        &IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder());

    std::vector<unsigned int> DOFs;

    int ndof = fe->GetDof();
    int dim = fe->GetDim();
    Array<int> vdofs;
    _pfes->GetElementDofs(elem_id, vdofs);

    DOFs.resize(vdofs.Size());
    for (int i = 0; i < vdofs.Size(); i++) {
      DOFs[i] = vdofs[i];
    }
    assert(DOFs.size() == 4);

    DenseMatrix temp_sigma = cond_att[eltrans->Attribute - 1];
    Matrix3D sigma_ele;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        sigma_ele(i, j) = temp_sigma(i, j);
      }
    }

    double elem_error = 0.;
    double elem_error_w = 0;
    double elem_error_e[solution.size() - 1];
    for (int i = 0; i < solution.size() - 1; i++)
      elem_error_e[i] = 0.;

    for (unsigned int q = 0; q < ir->GetNPoints(); q++) {

      const IntegrationPoint &ip = ir->IntPoint(q);
      eltrans->SetIntPoint(&ip);
      double WJ = eltrans->Weight() * ip.weight;

      DenseMatrix Dshape(ndof, 3);
      Dshape = 0.0;
      fe->CalcPhysDShape(*eltrans, Dshape);

      Point_ren dev_shape[4];
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
          dev_shape[i](j) = Dshape(i, j);
        }
      }

      Vector r_temp(3);
      r_temp = 0.0;
      eltrans->Transform(ip, r_temp);
      Point_ren r(r_temp[0], r_temp[1], r_temp[2]);

      DenseVector &W = solution[solution.size() - 1];
      Matrix43D &Wa = a[solution.size() - 1];

      Matrix13D grad_w_P;
      this->compute_P_value(r, Wa, grad_w_P);
      Matrix13D grad_w_h;
      Point_ren temp_grad_w_h;
      for (unsigned int t = 0; t < 4; t++)
        temp_grad_w_h = temp_grad_w_h + dev_shape[t] * W(DOFs[t]);
      for (int t = 0; t < 3; t++)
        grad_w_h(t) = temp_grad_w_h(t);

      elem_error_w += ((grad_w_P - grad_w_h) * sigma_ele *
                       (grad_w_P - grad_w_h).transpose())
                          .determinant() *
                      WJ;

      for (int i = 0; i < solution.size() - 1; i++) {
        DenseVector &Ui = solution[i];
        Matrix43D &Uia = a[i];
        Matrix13D grad_e_P;
        this->compute_P_value(r, Uia, grad_e_P);
        Matrix13D grad_e_h;
        Point_ren temp_grad_e_h;
        for (unsigned int t = 0; t < 4; t++)
          temp_grad_e_h = temp_grad_e_h + dev_shape[t] * Ui(DOFs[t]);
        for (int t = 0; t < 3; t++)
          grad_e_h(t) = temp_grad_e_h(t);

        elem_error_e[i] += ((grad_e_P - grad_e_h) * sigma_ele *
                            (grad_e_P - grad_e_h).transpose())
                               .determinant() *
                           WJ;
      }
    }

    for (int i = 0; i < solution.size() - 1; i++) {
      elem_error += std::sqrt(elem_error_w) * std::sqrt(elem_error_e[i]);
    }

    energy_norm[elem_id] = elem_error;
  }
  return;
}
