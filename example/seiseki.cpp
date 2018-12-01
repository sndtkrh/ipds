#include <iostream>
#include <ipds/ipds.hpp>

int main(){
  std::vector<std::vector<double>> seiseki(9);
  for(int i = 0; i < 166; i++) {
    for(int j = 0; j < 9; j++) {
      double p;
      std::cin >> p;
      seiseki[j].push_back(p);
    }
  }
  ipds::SymmetricMatrix seiseki_corel = ipds::corel_matrix(seiseki);
  auto eigen = seiseki_corel.get_eigen();
  double sum_of_lambda = 0;
  for(const auto & e : eigen) {
    sum_of_lambda += std::get<0>(e);
  }
  double cumulative_proportion = 0;
  for(const auto & e : eigen) {
    const auto & e_lambda = std::get<0>(e);
    const auto & e_vector = std::get<1>(e);
    cumulative_proportion += e_lambda / sum_of_lambda;
    /*
    std::cout
    << "Standard deviation: " << std::sqrt(e_lambda)
    << ", Proportion of Variance: " << e_lambda / sum_of_lambda
    << ", Cumulative Proportion: " << cumulative_proportion
    << std::endl;
    std::cout << "Eigen vector: ";
    for(int i = 0; i < e_vector.size(); i++) {
      std::cout << e_vector[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
    */
  }

  // svg test
  ipds::SVGcanvas svg(1000,1000);
  for(int i = 0; i < seiseki[0].size(); i++){
    svg.circles.emplace_back(seiseki[0][i] * 10, seiseki[1][i] * 10, 2);
  }
  svg.save("a.svg");
}
