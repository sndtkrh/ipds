#ifndef IPDS_SVG
#define IPDS_SVG

namespace ipds {

  struct Circle {
    double cx, cy, r;
    std::string color;
    Circle(double cx, double cy, double r, std::string color);
    std::string get_svg_string() const;
  };

  struct SVGcanvas {
    std::vector<Circle> circles;
    SVGcanvas(double w, double h);
    std::string get_svg_string() const;
    void save(std::string filename) const;
  private:
    double width, height;
  };
}
#endif
