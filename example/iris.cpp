#include <iostream>
#include <ipds/ipds.hpp>

int main(){
  int data_dim = 4;
  int data_n = 150;
  std::vector<std::vector<double>> iris(data_dim);
  for(std::size_t i = 0; i < data_n; i++) {
    for(std::size_t j = 0; j < iris.size(); j++) {
      double p;
      std::cin >> p;
      iris[j].push_back(p);
    }
  }

  for(std::size_t j = 0; j < iris.size(); j++) {
    iris[j] = ipds::standardize(iris[j]);
  }

  // Do PCA
  auto eigen = ipds::pca(iris);
  ipds::print_pca_result(eigen);

  // take 2-dim subspace and project data
  auto data_points = ipds::transpose(iris);
  auto basis = {std::get<1>(eigen[0]), std::get<1>(eigen[1])};
  int svgwh = 400;
  double scale = 50;
  ipds::SVGcanvas svg(svgwh,svgwh);
  int i = 0;
  for(const auto & p : data_points) {
    auto projected_p = ipds::project(p, basis);
    std::string color = (i < 50) ? "red" : ((i < 100) ? "green" : "blue");
    ipds::plot_point_2D({projected_p[0], projected_p[1]}, svg, scale, color);
    std::cout << "i=" << i++ << " (" << projected_p[0] << ", " << projected_p[1] << ") " << std::endl;
  }
  for(int i = 0; i < data_dim; i++) {
    std::vector<double> e(data_dim, 0);
    e[i] = 1;
    auto p = ipds::project(e, basis);
    ipds::plot_line_2D({0, 0}, {p[0], p[1]}, svg, scale);
  }
  svg.save("a.svg");
}
