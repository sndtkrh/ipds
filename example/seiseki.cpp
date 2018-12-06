#include <iostream>
#include <ipds/ipds.hpp>

int main(){
  std::size_t data_dim = 9;
  std::vector<std::vector<double>> seiseki(data_dim);
  for(std::size_t i = 0; i < 166; i++) {
    for(std::size_t j = 0; j < seiseki.size(); j++) {
      double p;
      std::cin >> p;
      seiseki[j].push_back(p);
    }
  }

  seiseki = ipds::standardize(ipds::polynomial_features(ipds::standardize(seiseki)));

  // Do PCA
  auto eigen = ipds::pca(seiseki);
  ipds::print_pca_result(eigen);

  // take 2-dim subspace and project data
  auto data_points = ipds::transpose(seiseki);
  auto basis = {std::get<1>(eigen[0]), std::get<1>(eigen[1])};
  int svgwh = 1000;
  double scale = 50;
  ipds::SVGcanvas svg(svgwh,svgwh);
  int i = 0;
  for(const auto & p : data_points) {
    auto projected_p = ipds::project(p, basis);
    ipds::plot_point_2D({projected_p[0], projected_p[1]}, svg, scale, "black", 2);
    std::cout << "i=" << i++ << " (" << projected_p[0] << ", " << projected_p[1] << ") " << std::endl;
  }
  for(int i = 0; i < seiseki.size(); i++) {
    std::vector<double> e(seiseki.size(), 0);
    e[i] = 1;
    auto p = ipds::project(e, basis);
    ipds::plot_line_2D({0, 0}, {p[0], p[1]}, svg, scale * 10, "blue");
    svg.texts.emplace_back(std::to_string(i), svgwh / 2 + p[0] * scale * 10, svgwh / 2 + p[1] * scale * 10, 5, "blue");
  }
  svg.save("seiseki.svg");
}
