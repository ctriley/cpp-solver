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
  private:
    static void updateInverse(const Eigen::VectorXd &alpha, long leaving_index,
        Eigen::MatrixXd &Abeta_inverse);
    static long findLeavingValue(const Eigen::VectorXd &x_beta,
        const Eigen::VectorXd &alpha);
    [[nodiscard]] static std::pair<long, double> findEnteringValue(OBJ_FUNC mode,
        const Eigen::VectorXd &rc);
    template <typename T>
    std::pair<bool, int > findInVector(const std::vector<T>  & vecOfElements,
        const T  & element);
    template <typename T>
    void printVector(const std::vector<T> &vec);

  public:
    Simplex(OBJ_FUNC mode, const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
        const Eigen::VectorXd &c, std::vector<long> &basis);
  };
}

#endif //CPP_SOLVER_SIMPLEX_HPP
