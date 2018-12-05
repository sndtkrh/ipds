#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <tuple>
#include "svg/svg.hpp"

namespace ipds {
  Circle::Circle(double cx, double cy, double r, std::string color)
    : cx(cx), cy(cy), r(r), color(color) {}
  std::string Circle::get_svg_string() const {
    std::stringstream ss;
    ss << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << r
      << "\" fill=\"" << color << "\" ></circle>";
    return ss.str();
  }

  Line::Line(double x1, double y1, double x2, double y2, std::string color, double stroke_width)
    : x1(x1), y1(y1), x2(x2), y2(y2), color(color), stroke_width(stroke_width) {}
  std::string Line::get_svg_string() const {
    std::stringstream ss;
    ss << "<line x1=\"" << x1 << "\" y1=\"" << y1
      << "\" x2=\"" << x2 << "\" y2=\"" << y2
      << "\" stroke=\"" << color
      << "\" stroke-width=\"" << stroke_width << "\" ></line>";
    return ss.str();
  }

  Text::Text(std::string text, double x, double y, double font_size, std::string color)
    : text(text), x(x), y(y), font_size(font_size), color(color) {}
  std::string Text::get_svg_string() const {
    std::stringstream ss;
    ss << "<text x=\"" << x << "\" y=\"" << y
      << "\" font-size=\"" << font_size
      << "\" fill=\"" << color << "\" >"
      << text
      << "</text>";
    return ss.str();
  }

  SVGcanvas::SVGcanvas(double w, double h) : width(w), height(h) {}
  std::string SVGcanvas::get_svg_string() const {
    std::stringstream ss;
    ss << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
    << "\" height=\"" << height
    << "\" viewBox=\"0 0 " << width << " " << height << "\">\n";

    for(const Line & line : lines) {
      ss << line.get_svg_string() << "\n";
    }
    for(const Circle & circle : circles) {
      ss << circle.get_svg_string() << "\n";
    }
    for(const Text & text : texts) {
      ss << text.get_svg_string() << "\n";
    }

    ss << "</svg>\n";
    return ss.str();
  }

  void SVGcanvas::save(std::string filename) const {
    std::ofstream ofs(filename);
    ofs << get_svg_string();
  }

  std::tuple<double, double> SVGcanvas::get_canvas_size() const {
    return {width, height};
  }
}
