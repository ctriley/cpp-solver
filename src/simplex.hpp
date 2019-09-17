//
// Created by connor on 9/5/19.
//

#ifndef CPP_SOLVER_SIMPLEX_HPP
#define CPP_SOLVER_SIMPLEX_HPP

#include <vector>
#include <Eigen/Dense>

namespace cppsolver {
  typedef enum {MAXIMIZE, MINIMIZE} OBJ_FUNC;
  typedef enum {PRIMAL, DUAL, PHASE_I} SIMPLEX_MODE;

  class Simplex {
  private:
    void runPrimal(OBJ_FUNC mode, const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
            const Eigen::VectorXd &c, std::vector<long> &basis);
    void runDual(OBJ_FUNC mode, const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
            const Eigen::VectorXd &c, std::vector<long> &basis);
    void runPhaseI(OBJ_FUNC mode, Eigen::MatrixXd &A, Eigen::VectorXd &b,
            Eigen::VectorXd &c, std::vector<long> &basis);
    SIMPLEX_MODE chooseSimplexMode(const Eigen::MatrixXd &A,
        const Eigen::VectorXd &b, const Eigen::VectorXd &c, std::vector<long> &basis);
    static void updateInverse(const Eigen::VectorXd &alpha, long leaving_index,
        Eigen::MatrixXd &Abeta_inverse);
    static long findLeavingValue(const Eigen::VectorXd &x_beta,
        const Eigen::VectorXd &alpha, int multiplier=1);
    [[nodiscard]] static std::pair<long, double> findEnteringValue(OBJ_FUNC mode,
        const Eigen::VectorXd &rc);
    void createNonBasis(const Eigen::MatrixXd &A,
        const std::vector<long> &basis, std::vector<long> &non_basis);
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
