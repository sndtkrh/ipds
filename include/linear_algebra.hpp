#ifndef IPDS_MATRIX
#define IPDS_MATRIX

namespace ipds {

  class Matrix {
  public:
    Matrix(const std::vector<std::vector<double>> & mat_);
    const std::vector<double> & operator[] (int i) const;
    std::vector<double> operator() (const std::vector<double> & v) const;
    std::size_t rows() const;
    std::size_t cols() const;
  protected:
    std::vector<std::vector<double>> mat;
    std::size_t row, col;
  };

  class SquareMatrix : public Matrix {
  public:
    SquareMatrix(const std::vector<std::vector<double>> & mat_);
  };

  class SymmetricMatrix : public SquareMatrix {
  public:
    SymmetricMatrix(const std::vector<std::vector<double>> & mat_);
    std::vector<std::tuple<double, std::vector<double>>> get_eigen() const;
  };

  double inner_product(const std::vector<double> & v, const std::vector<double> & u);
  double norm(const std::vector<double> & v);
  std::vector<std::vector<double>> orthonormalization(const std::vector<std::vector<double>> & v);
}

#endif
