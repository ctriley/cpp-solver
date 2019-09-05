#include <iostream>

#include <Eigen/Dense>
#include <random>

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"

#include "config.hpp"
#include "simplex.hpp"

int main(int argc, char* argv[]) {

  std::string app_name = boost::filesystem::basename(argv[0]);
  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()
      ("help,h", "produce a help message")
      ("version,v", "output the version number");
  po::positional_options_description pd;
  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(desc).positional(
        pd).run(), vm);
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      return 0;
    }
    if (vm.count("version")) {
      std::cout << PROJECT_VERSION << std::endl;
      return 0;
    }
    po::notify(vm);

  } catch(po::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    std::cerr << desc << std::endl;
    return 1;
  }
  std::cout << "Hello, World!" << std::endl;



  // if just swapping two rols or two columns you can do this:
  Eigen::MatrixXi M = Eigen::MatrixXi::Random(3,3);
  std::cout << M << std::endl;
  M.row(1).swap(M.row(2));
  std::cout << "**************************" << std::endl;
  std::cout << M << std::endl;
  std::cout << "**************************" << std::endl;

  // a permutation matrix (row, col) = (origin, destination)
  // if you do A * perm you permute columns
  // perm * A permutes the rows

  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(3);
  perm.setIdentity();
  Eigen::MatrixXi M_perm = M * perm; // permute columns
  std::cout << M << std::endl;
  std::cout << "**************************" << std::endl;
  std:: cout << M_perm << std::endl;

  M_perm = perm * M; // permute rows


  // intializing values
  Eigen::VectorXd c(5);
  c << -2, 7, 14, 2, -4;
  Eigen::VectorXd b(2);
  b << 0, 7;
  Eigen::MatrixXd A(2, 5);
  A << -2, -2, 3, 1, -2,
  1, 4, 1, -1, -1;

  std::vector<int> primal_basis = {0, 2};
  std::vector<int> dual_basis = {3, 4};
  cppsolver::Simplex(cppsolver::OBJ_FUNC::MINIMIZE, A, b, c, primal_basis);
  return 0;
}