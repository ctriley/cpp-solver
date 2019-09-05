//
// Created by connor on 9/5/19.
//

#include "simplex.hpp"

namespace cppsolver {

  Simplex::Simplex(OBJ_FUNC mode, const Eigen::MatrixXd &A,
      const Eigen::VectorXd &b, const Eigen::VectorXd &c, std::vector<int> &basis) {
    Eigen::MatrixXd Abeta_init = A(Eigen::all, basis);
    Eigen::VectorXd cbeta = c(Eigen::all, basis);
    Eigen::MatrixXd Abeta_inverse = Abeta_init.inverse();

    // compute values needed every single time
    Eigen::VectorXd x_beta = Abeta_inverse * b;
    Eigen::VectorXd rc = cbeta * Abeta_inverse;

    // find entering value
    long worst_rc_index = 0;
    long worst_rc = 0;
    for(long i = 0; i < rc.size(); ++i) {
      if(rc[i] < worst_rc) {
        worst_rc = rc[i];
        worst_rc_index = i;
      }
    }
    if(worst_rc > 0) {
      // stop
    }

    // pick leaving value

  }


  /*template<typename Derived>
  void Simplex::primal(int mode, const Eigen::VectorXd &c,
      const Eigen::MatrixBase<Derived> &A, const Eigen::VectorXd &b,
      const std::vector<int> beta) {
    Eigen::PermutationMatrix<Eigen::Dynamic> permutation_matrix;
    permutation_matrix.setIdentity();
    permutation_matrix.
    Eigen::MatrixX2d cbeta = c
    Eigen::MatrixX2d ybar = c[b];
  }*/
}