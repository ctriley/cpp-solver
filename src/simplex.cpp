//
// Created by connor on 9/5/19.
//

#include "simplex.hpp"

#include <iostream>
#include <algorithm>

namespace cppsolver {

  Simplex::Simplex(OBJ_FUNC mode, const Eigen::MatrixXd &A,
      const Eigen::VectorXd &b, const Eigen::VectorXd &c, std::vector<long> &basis) {
    int iteration = 0;
    auto simplex_mode = chooseSimplexMode(A, b, c, basis);
    if (simplex_mode == SIMPLEX_MODE::PRIMAL) {
      runPrimal(mode, A, b, c, basis);
    } else if (simplex_mode == SIMPLEX_MODE::DUAL) {
      runDual(mode, A, b, c, basis);
    } else if (simplex_mode == SIMPLEX_MODE::PHASE_I) {
      Eigen::MatrixXd A_copy(A);
      Eigen::VectorXd b_copy(b);
      Eigen::VectorXd c_copy(c);
      runPhaseI(mode, A_copy, b_copy, c_copy, basis);
      runPrimal(mode, A_copy, b_copy, c_copy, basis);
    }
  }


  void Simplex::runDual(OBJ_FUNC mode, const Eigen::MatrixXd &A,
    const Eigen::VectorXd &b, const Eigen::VectorXd &c, std::vector<long> &basis) {

    std::vector<long> non_basis;
    createNonBasis(A, basis, non_basis);
    int iteration = 0;
    Eigen::MatrixXd Abeta = A(Eigen::all, basis);

    Eigen::MatrixXd Abeta_inverse = Abeta.inverse();
    while(true) {
      // create permutation matrix to permute A and c
      std::vector<long> index_vector = basis;
      index_vector.insert(index_vector.end(), non_basis.begin(), non_basis.end());
      Eigen::VectorXi indices(A.cols());
      for(long i = 0; i < indices.size(); ++i) {
        indices[i] = index_vector[i];
      }
      Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
      perm.indices() = indices;
      Abeta = A(Eigen::all, basis);
      Eigen::MatrixXd Abeta_inverse_test = Abeta.inverse();
      std::cout << Abeta_inverse << std::endl;
      std::cout << "*******************" << std::endl;
      std::cout << Abeta_inverse_test << std::endl;
      Eigen::VectorXd c_beta = c(basis, Eigen::all);
      Eigen::VectorXd pi_beta = c_beta.transpose() * Abeta_inverse;
      Eigen::MatrixXd A_eta = A(Eigen::all, non_basis);
      Eigen::VectorXd c_eta = c(non_basis, Eigen::all);
      Eigen::MatrixXd A_perm = A * perm;
      Eigen::MatrixXd c_perm = c.transpose() * perm;
      Eigen::VectorXd pi = Eigen::VectorXd::Zero(c_perm.size());
      pi.topRows(basis.size()) = pi_beta.topRows(basis.size());
      Eigen::VectorXd temp = (c.transpose() * perm);
      Eigen::VectorXd temp2 = (pi.transpose() * A_perm.transpose());
      Eigen::VectorXd d = (c.transpose() * perm);
      d.topRows(basis.size()) = d.topRows(basis.size()) - temp2;
      Eigen::VectorXd x_beta = Abeta_inverse * b;
      if(x_beta.minCoeff() >= 0) {
        std::cout << "Num iterations: " << iteration << std::endl;
        std::cout << "Objective: " << b.transpose() * pi_beta << std::endl;
        std::cout << "X values: " << pi_beta << std::endl;
        break;
      }
      // now choose leaving and entering values
      long leaving_index = 0;
      double lowest_value = std::numeric_limits<double>::max();
      for(long i = 0; i < x_beta.size(); ++i) {
        if(x_beta[i] < lowest_value) {
          leaving_index = i;
          lowest_value = x_beta[i];
        }
      }
      Eigen::VectorXd t = pi_beta.transpose() * A_eta - c_eta.transpose();
      Eigen::VectorXd v = Abeta_inverse.row(leaving_index) * A_eta;
      if(v.minCoeff() >= 0) {
        std::cout << "Infeasible" << std::endl;
        break;
      }
      long entering_index = 0;
      double min_ratio = std::numeric_limits<double>::max();
      for(long i = 0; i < non_basis.size(); ++i) {
        if(v[i] < 0) {
          if(t[i] / v[i] < min_ratio) {
            min_ratio = t[i] / v[i];
            entering_index = i;
          }
        }
      }
      Eigen::VectorXd d_eta = d.bottomRows(non_basis.size());
      entering_index = non_basis.at(entering_index);
      Eigen::VectorXd y = A.row(entering_index);
      updateInverse(v, leaving_index, Abeta_inverse);
      leaving_index = basis.at(leaving_index);
      // update leaving index to original A index
      auto p1 = findInVector<long>(basis, leaving_index);
      auto p2 = findInVector<long>(non_basis, entering_index);
      basis[p1.second] = entering_index;
      non_basis[p2.second] = leaving_index;
      std::cout << "basis: ";
      printVector(basis);
      std::cout << "nonbasis: ";
      printVector(non_basis);
      ++iteration;
    }
  }

  void Simplex::runPhaseI(OBJ_FUNC mode, Eigen::MatrixXd &A,
    Eigen::VectorXd &b, Eigen::VectorXd &c, std::vector<long> &basis) {

  }

  void Simplex::runPrimal(OBJ_FUNC mode, const Eigen::MatrixXd &A,
    const Eigen::VectorXd &b, const Eigen::VectorXd &c, std::vector<long> &basis) {
    std::vector<long> non_basis;
    createNonBasis(A, basis, non_basis);
    Eigen::MatrixXd Abeta = A(Eigen::all, basis);
    Eigen::MatrixXd Abeta_inverse = Abeta.inverse();
    int iteration = 0;
    while(true) {
      // get slices
      Eigen::MatrixXd A_eta = A(Eigen::all, non_basis);
      Eigen::VectorXd cbeta = c(basis, Eigen::all);
      Eigen::VectorXd c_eta = c(non_basis, Eigen::all);
      // compute values
      Eigen::VectorXd x_beta = Abeta_inverse * b;
      Eigen::VectorXd ybar = cbeta.transpose() * Abeta_inverse;
      Eigen::VectorXd rc = c_eta - (A_eta.transpose() * ybar);
      // find entering value
      auto [entering_index, worst_rc] = findEnteringValue(mode, rc);
      if (worst_rc > 0) {
        std::cout << "Num iterations: " << iteration << std::endl;
        std::cout << "Objective: " << cbeta.transpose() * x_beta << std::endl;
        std::cout << "X values: " << x_beta << std::endl;
        printVector(basis);
        break;
        // stop
      }
      // transform entering index to original A index
      entering_index = non_basis.at(entering_index);
      // find leaving index
      Eigen::VectorXd entering_vector = A.col(entering_index);
      Eigen::VectorXd alpha = Abeta_inverse * entering_vector;
      auto leaving_index = findLeavingValue(x_beta, alpha);
      if (leaving_index < 0) {
        std::cout << "problem is unbounded" << std::endl;
        break;
      }
      // update inverse
      updateInverse(alpha, leaving_index, Abeta_inverse);
      // update leaving index to original A index
      leaving_index = basis.at(leaving_index);
      auto p1 = findInVector<long>(basis, leaving_index);
      auto p2 = findInVector<long>(non_basis, entering_index);
      basis[p1.second] = entering_index;
      non_basis[p2.second] = leaving_index;
      ++iteration;
    }
  }

  SIMPLEX_MODE Simplex::chooseSimplexMode(const Eigen::MatrixXd &A,
      const Eigen::VectorXd &b, const Eigen::VectorXd &c, std::vector<long> &basis) {
    std::vector<int> non_basis;
    for (long i = 0; i < A.cols(); ++i) {
      if (!findInVector(basis, i).first) {
        non_basis.push_back(i);
      }
    }

    // create permutation matrix to permute A and c
    std::vector<long> index_vector = basis;
    index_vector.insert(index_vector.end(), non_basis.begin(), non_basis.end());
    Eigen::VectorXi indices(A.cols());
    for(long i = 0; i < indices.size(); ++i) {
      indices[i] = index_vector[i];
    }
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
    perm.indices() = indices;

    // compute primal and dual variable values
    Eigen::MatrixXd A_beta_inverse = A(Eigen::all, basis).inverse();
    Eigen::VectorXd y_beta = c(basis, Eigen::all).transpose() * A_beta_inverse;
    Eigen::VectorXd x_beta = A_beta_inverse * b;

    // x_beta and y_beta are the primal and dual values
    if(x_beta.minCoeff() >= 0) {
      // primal simplex
      return SIMPLEX_MODE::PRIMAL;
    } else if(((c.transpose() * perm) - y_beta.transpose() * (A * perm)).minCoeff() >= 0) {
      // dual simplex
      return SIMPLEX_MODE::DUAL;
    }
    return SIMPLEX_MODE::PHASE_I;

  }

  void Simplex::updateInverse(const Eigen::VectorXd &alpha, long leaving_index,
      Eigen::MatrixXd &Abeta_inverse) {
    Eigen::VectorXd temp(alpha.size());
    for(long i = 0; i < alpha.size(); ++i) {
      temp[i] = - alpha[i] / alpha[leaving_index];
    }
    temp[leaving_index] = 1 / alpha[leaving_index];
    Eigen::MatrixXd E = Eigen::MatrixXd::Identity(alpha.size(), Abeta_inverse.rows());
    E.col(leaving_index) = temp;
    Abeta_inverse = E * Abeta_inverse;
  }


  long Simplex::findLeavingValue(const Eigen::VectorXd &x_beta,
      const Eigen::VectorXd &alpha, int multiplier) {
    // pick leaving value
    long leaving_index = -1;
    double lowest_theta = std::numeric_limits<double>::max();
    for (long i = 0; i < x_beta.size(); ++i) {
      if (alpha[i] > 0) {
        double theta = multiplier * x_beta[i] / alpha[i];
        if (theta < lowest_theta) {
          lowest_theta = theta;
          leaving_index = i;
        }
      }
    }
    return leaving_index;
  }

  std::pair<long, double> Simplex::findEnteringValue(OBJ_FUNC mode,
      const Eigen::VectorXd &rc) {
    long entering_index = 0;
    double worst_rc = rc[0];
    for (long i = 0; i < rc.size(); ++i) {
      if (mode == OBJ_FUNC::MINIMIZE) {
        if (rc[i] < worst_rc) {
          worst_rc = rc[i];
          entering_index = i;
        }
      } else {
        if (rc[i] > worst_rc) {
          worst_rc = rc[i];
          entering_index = i;
        }
      }
    }
    return {entering_index, worst_rc};
  }

  template <typename T>
  std::pair<bool, int > Simplex::findInVector(const std::vector<T>  &vecOfElements,
      const T  &element) {
    std::pair<bool, int> result;
    // Find given element in vector
    auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);
    if (it != vecOfElements.end()) {
      result.second = distance(vecOfElements.begin(), it);
      result.first = true;
    } else {
      result.first = false;
      result.second = -1;
    }
    return result;
  }

  void Simplex::createNonBasis(const Eigen::MatrixXd &A,
      const std::vector<long> &basis, std::vector<long> &non_basis) {
    for (long i = 0; i < A.cols(); ++i) {
      if (!findInVector(basis, i).first) {
        non_basis.push_back(i);
      }
    }
  }

  template <typename T>
  void Simplex::printVector(const std::vector<T> &vec) {
    for(auto i : vec) {
      std::cout << i << ", ";
    }
    std::cout << std::endl;
  }
}