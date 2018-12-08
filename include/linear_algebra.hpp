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

    /*
    * get_eigen:
    *   return value: {(e_val, e_vec)}
    *     - e_val is the eigen max_value.
    *     - e_vec is the eigen vector.
    *     - the returned vector is sorted: |e_val_0| >= ... >= |e_val_n|.
    */
    std::vector<std::tuple<double, std::vector<double>>> get_eigen_yacobi() const;
    std::vector<std::tuple<double, std::vector<double>>> get_eigen_householder_binarysearch() const;

    /*
    * tridiagonalize:
    *   return value: (T, {v})
    *     - T is the tridiagonal matrix of mat.
    *     - vs are the vector representing householder transformations.
    */
    std::tuple<SymmetricMatrix, std::vector<std::vector<double>>> tridiagonalize() const;
  };

  /*
  * householder:
  *   argument: u, idx
  *     - u is a vector that you want to transform.
  *     - idx is for tridiagonalize.
  *   return value: (s, v)
  *     - s stisfies [s 0 ... 0]^T = Qu whre Q := I - 2v(v^T).
  *     - v is the normal vector of householder transformation.
  */
  std::tuple<double, std::vector<double>> householder(const std::vector<double> & u, const std::size_t idx = 0);

  std::vector<double> operator + (const std::vector<double> & a, const std::vector<double> & b);
  std::vector<double> operator * (const double lambda, const std::vector<double> & a);
  std::vector<double> operator * (const Matrix & mat, const std::vector<double> & a);

  double inner_product(const std::vector<double> & v, const std::vector<double> & u);
  double norm(const std::vector<double> & v);
  std::vector<std::vector<double>> orthonormalize(const std::vector<std::vector<double>> & v);
  std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> & a);
  std::vector<double> project(const std::vector<double> & point, const std::vector<std::vector<double>> & basis_of_subspace);
  std::vector<std::vector<double>> project(const std::vector<std::vector<double>> & points, const std::vector<std::vector<double>> & basis_of_subspace);
}

#endif
