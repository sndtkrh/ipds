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

  struct Text {
    std::string text;
    double x, y;
    double font_size;
    std::string color;
    Text(std::string text, double x, double y, double font_size, std::string color = "black");
    std::string get_svg_string() const;
  };

  struct SVGcanvas {
    std::vector<Line> lines;
    std::vector<Circle> circles;
    std::vector<Text> texts;
    SVGcanvas(double w, double h);
    std::string get_svg_string() const;
    void save(std::string filename) const;
    std::tuple<double, double> get_canvas_size() const;
  private:
    double width, height;
  };

  void plot_point_2D (std::tuple<double, double> point, SVGcanvas & svg, double scale = 1, std::string color = "black", double point_radius = 1);
  void plot_line_2D (std::tuple<double, double> point_begin, std::tuple<double, double> point_end, SVGcanvas & svg, double scale = 1, std::string color = "black", double stroke_width = 1);
}
#endif
