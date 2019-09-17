//
// Created by connor on 9/5/19.
//

#include "simplex.hpp"
#include "gtest/gtest.h"


namespace cppsolver {
  class SimplexTest : public ::testing::Test {
  protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    SimplexTest() {
      // You can do set-up work for each test here.
    }

    ~SimplexTest() override {
      // You can do clean-up work that doesn't throw exceptions here.
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    void SetUp() override {
      // Code here will be called immediately after the constructor (right
      // before each test).
    }

    void TearDown() override {
      // Code here will be called immediately after each test (right
      // before the destructor).
    }

    // Objects declared here can be used by all tests in the test suite for Foo.
  };

  // Tests that the Foo::Bar() method does Abc.
  TEST_F(SimplexTest, PrimalTest) {

    // intializing values
    Eigen::VectorXd c(5);
    c << -2, 7, 14, 2, -4;
    Eigen::VectorXd b(2);
    b << 0, 7;
    Eigen::MatrixXd A(2, 5);
    A << -2, -2, 3, 1, -2,
        1, 4, 1, -1, -1;
    std::vector<long> primal_basis = {0, 2};

    Simplex simplex(OBJ_FUNC::MINIMIZE, A, b, c, primal_basis);
    ASSERT_DOUBLE_EQ(simplex.getObjective(), 24.5);
  }

  TEST_F(SimplexTest, DualTest) {
    // intializing values
    Eigen::VectorXd c(5);
    c << -2, 7, 14, 2, -4;
    Eigen::VectorXd b(2);
    b << 0, 7;
    Eigen::MatrixXd A(2, 5);
    A << -2, -2, 3, 1, -2,
        1, 4, 1, -1, -1;
    std::vector<long> dual_basis = {3, 4};
    Simplex simplex(OBJ_FUNC::MINIMIZE, A, b, c, dual_basis);
    ASSERT_DOUBLE_EQ(simplex.getObjective(), 24.5);

  }

  TEST_F(SimplexTest, PhaseITEST) {
    // intializing values
    Eigen::VectorXd c(5);
    c << -2, 7, 14, 2, -4;
    Eigen::VectorXd b(2);
    b << 0, 7;
    Eigen::MatrixXd A(2, 5);
    A << -2, -2, 3, 1, -2,
        1, 4, 1, -1, -1;
    std::vector<long> not_a_basis = {1, 4};
    cppsolver::Simplex simplex(cppsolver::OBJ_FUNC::MINIMIZE, A, b, c,
                       not_a_basis);
    ASSERT_DOUBLE_EQ(simplex.getObjective(), 24.5);

  }

  TEST_F(SimplexTest, PhaseIEmptyTest) {
    // intializing values
    Eigen::VectorXd c(5);
    c << -2, 7, 14, 2, -4;
    Eigen::VectorXd b(2);
    b << 0, 7;
    Eigen::MatrixXd A(2, 5);
    A << -2, -2, 3, 1, -2,
        1, 4, 1, -1, -1;
    std::vector<long> empty_basis = {};
    cppsolver::Simplex simplex(cppsolver::OBJ_FUNC::MINIMIZE, A, b, c,
        empty_basis);
    ASSERT_DOUBLE_EQ(simplex.getObjective(), 24.5);
  }
}
