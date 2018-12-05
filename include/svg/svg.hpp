#ifndef IPDS_SVG
#define IPDS_SVG

namespace ipds {

  struct Circle {
    double cx, cy, r;
    std::string color;
    Circle(double cx, double cy, double r, std::string color = "black");
    std::string get_svg_string() const;
  };

  struct Line {
    double x1, y1, x2, y2;
    std::string color;
    double stroke_width;
    Line(double x1, double y1, double x2, double y2, std::string color = "black", double stroke_width = 1);
    std::string get_svg_string() const;
  };

  struct SVGcanvas {
    std::vector<Circle> circles;
    std::vector<Line> lines;
    SVGcanvas(double w, double h);
    std::string get_svg_string() const;
    void save(std::string filename) const;
    std::tuple<double, double> get_canvas_size() const;
  private:
    double width, height;
  };
}
#endif
