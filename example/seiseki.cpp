#include <iostream>
#include <ipds/ipds.hpp>

int main(){
  std::vector<std::vector<double>> seiseki(9);
  for(std::size_t i = 0; i < 166; i++) {
    for(std::size_t j = 0; j < 9; j++) {
      double p;
      std::cin >> p;
      seiseki[j].push_back(p);
    }
  }

  for(std::size_t j = 0; j < seiseki.size(); j++) {
    seiseki[j] = ipds::standardize(seiseki[j]);
  }

  // Do PCA
  auto eigen = ipds::pca(seiseki);
  double sum_of_lambda = 0;
  for(const auto & e : eigen) {
    sum_of_lambda += std::get<0>(e);
  }
  double cumulative_proportion = 0;
  for(const auto & e : eigen) {
    const auto & e_lambda = std::get<0>(e);
    const auto & e_vector = std::get<1>(e);
    cumulative_proportion += e_lambda / sum_of_lambda;
    std::cout
    << "\x1b[36mStandard deviation\x1b[39m: " << std::sqrt(e_lambda)
    << ", \x1b[36mProportion of Variance\x1b[39m: " << e_lambda / sum_of_lambda
    << ", \x1b[36mCumulative Proportion\x1b[39m: " << cumulative_proportion
    << std::endl;
    std::cout << "\x1b[36mEigen vector\x1b[39m: ";
    for(std::size_t i = 0; i < e_vector.size(); i++) {
      std::cout << e_vector[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }

  // take 2-dim subspace and project data
  auto data_points = ipds::transpose(seiseki);
  auto basis = {std::get<1>(eigen[0]), std::get<1>(eigen[1])};
  int svgwh = 1000;
  ipds::SVGcanvas svg(svgwh,svgwh);
  int i = 0;
  for(const auto & p : data_points) {
    auto projected_p = ipds::project(p, basis);
    svg.circles.emplace_back(svgwh/2 + projected_p[0] * 50, svgwh/2 + projected_p[1] * 50, 2, ((i == 3) || (i == 129)) ? "green" : "black");
    std::cout << "i=" << i++ << " (" << projected_p[0] << ", " << projected_p[1] << ") " << std::endl;
  }
  for(int i = 0; i < 9; i++) {
    std::vector<double> e(9, 0);
    e[i] = 5;
    auto p = ipds::project(e, basis);
    svg.circles.emplace_back(svgwh/2 + p[0] * 50, svgwh/2 + p[1] * 50, 2, "red");
  }
  svg.circles.emplace_back(svgwh/2, svgwh/2, 2, "blue");
  svg.save("a.svg");
}
