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
  // intializing values
  Eigen::VectorXd c(5);
  c << -2, 7, 14, 2, -4;
  Eigen::VectorXd b(2);
  b << 0, 7;
  Eigen::MatrixXd A(2, 5);
  A << -2, -2, 3, 1, -2,
  1, 4, 1, -1, -1;

  std::vector<long> primal_basis = {0, 2};
  std::vector<long> dual_basis = {3, 4};
  std::vector<long> not_a_basis = {1, 4};
  cppsolver::Simplex(cppsolver::OBJ_FUNC::MINIMIZE, A, b, c,
      dual_basis);
  return 0;
}
