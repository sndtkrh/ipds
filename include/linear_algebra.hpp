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

  std::vector<double> operator + (const std::vector<double> & a, const std::vector<double> & b);
  std::vector<double> operator * (const double lambda, const std::vector<double> & a);

  double inner_product(const std::vector<double> & v, const std::vector<double> & u);
  double norm(const std::vector<double> & v);
  std::vector<std::vector<double>> orthonormalize(const std::vector<std::vector<double>> & v);
  std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> & a);
  std::vector<double> project(const std::vector<double> & point, const std::vector<std::vector<double>> & basis_of_subspace);
  std::vector<std::vector<double>> project(const std::vector<std::vector<double>> & points, const std::vector<std::vector<double>> & basis_of_subspace);
}

#endif
