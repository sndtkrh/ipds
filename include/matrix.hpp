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

  double covar(const std::vector<double> & x, const std::vector<double> & y);
  SymmetricMatrix var_covar_matrix(const std::vector<std::vector<double>> & data);
  SymmetricMatrix corel_matrix(const std::vector<std::vector<double>> & data);
}

#endif
