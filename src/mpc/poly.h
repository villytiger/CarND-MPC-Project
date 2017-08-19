#ifndef MPC_POLY_H_
#define MPC_POLY_H_

#include <cppad/cppad.hpp>

#include "Eigen/Core"
#include "Eigen/QR"


template<typename T>
class Poly {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;
  
  static Poly Fit(const Vector& xvals, const Vector& yvals, int order) {
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);

    Eigen::MatrixXd A(xvals.size(), order + 1);
    A.col(0).array() = 1;
    for (int i = 0; i != order; ++i)
      A.col(i + 1) = A.col(i).array() * xvals.array();

    return Poly(A.householderQr().solve(yvals));
  }

  Poly(const Vector& coeffs)
    : coeffs_(coeffs) {}


  T operator()(T x, unsigned k = 0) const {
    return CppAD::Poly(k, coeffs_, x);
  }

  const Vector& coeffs() const {
    return coeffs_;
  }
  
private:  
  Vector coeffs_;
};

#endif // MPC_POLY_H_
