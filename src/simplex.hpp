//
// Created by connor on 9/5/19.
//

#ifndef CPP_SOLVER_SIMPLEX_HPP
#define CPP_SOLVER_SIMPLEX_HPP

#include <vector>

#include <Eigen/Dense>

namespace cppsolver {
  typedef enum {MAXIMIZE, MINIMIZE} OBJ_FUNC;

  class Simplex {
  public:
    Simplex(OBJ_FUNC mode, const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
        const Eigen::VectorXd &c, std::vector<int> &basis);
  };
}

#endif //CPP_SOLVER_SIMPLEX_HPP
