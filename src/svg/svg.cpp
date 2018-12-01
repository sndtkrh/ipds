#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include "svg/svg.hpp"

namespace ipds {
  Circle::Circle(double cx, double cy, double r) : cx(cx), cy(cy), r(r) {}
  std::string Circle::get_svg_string() const {
    std::stringstream ss;
    ss << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << r << "\" ></circle>";
    return ss.str();
  }

  SVGcanvas::SVGcanvas(double w, double h) : width(w), height(h) {}

  std::string SVGcanvas::get_svg_string() const {
    std::stringstream ss;
    ss << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
    << "\" height=\"" << height
    << "\" viewBox=\"0 0 " << width << " " << height << "\">\n";

    for(const Circle & circle : circles) {
      ss << circle.get_svg_string() << "\n";
    }

    ss << "</svg>\n";
    return ss.str();
  }

  void SVGcanvas::save(std::string filename) const {
    std::ofstream ofs(filename);
    ofs << get_svg_string();
  }

}
