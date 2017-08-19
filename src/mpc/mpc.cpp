#include "mpc.h"

#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>

#include "Eigen/Core"

#include "poly.h"

using std::vector;
using Eigen::VectorXd;
using CppAD::AD;

constexpr const size_t kStateSize = 6;
constexpr const size_t kExStateSize = 8;

constexpr const size_t kSteps = 10;
constexpr const double kTimeDelta = 0.15;

constexpr const size_t kIndexX = 0;
constexpr const size_t kIndexY = 1;
constexpr const size_t kIndexPsi = 2;
constexpr const size_t kIndexV = 3;
constexpr const size_t kIndexCte = 4;
constexpr const size_t kIndexEpsi = 5;
constexpr const size_t kIndexDelta = 6;
constexpr const size_t kIndexA = 7;

constexpr const double kCostCte = 1;
constexpr const double kCostEpsi = 100;
constexpr const double kCostV = 1;
constexpr const double kCostDelta = 1;
constexpr const double kCostA = 1;
constexpr const double kCostDiffDelta = 5;
constexpr const double kCostDiffA = 1;

constexpr const double kLf = 2.67;

template<typename To, typename From>
To convert(const From& from) {
  To to(from.size());
  for (decltype(from.size()) i = 0; i != from.size(); ++i)
    to[i] = from[i];
  return to;
}

class FgEval {
public:
  typedef Eigen::Matrix<CppAD::AD<double>, Eigen::Dynamic, 1> ADvector;
  typedef Eigen::Matrix<CppAD::AD<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
  
  FgEval(const Poly<double>& poly, double speed)
    : poly_(poly.coeffs().cast<CppAD::AD<double>>())
    , speed_(speed) {
  }

  void operator()(ADvector& fg_vector, const ADvector& vars_vector) {
    Eigen::Map<Matrix> fg(fg_vector.data() + 1 + kExStateSize, kSteps - 1, kStateSize);
    Eigen::Map<const Matrix> vars(vars_vector.data(), kSteps, kExStateSize);

    fg_vector[0] = 0;

    fg_vector[0] += kCostCte * vars.col(kIndexCte).squaredNorm();
    fg_vector[0] += kCostEpsi * vars.col(kIndexEpsi).squaredNorm();
    fg_vector[0] += kCostV * (vars.col(kIndexV) - ADvector::Constant(kSteps, speed_)).squaredNorm();

    fg_vector[0] += kCostDelta * vars.bottomRows(kSteps - 1).col(kIndexDelta).squaredNorm();
    fg_vector[0] += kCostA * vars.bottomRows(kSteps - 1).col(kIndexA).squaredNorm();

    for (size_t i = 1; i < kSteps; ++i) {
      fg_vector[0] += kCostDiffDelta * CppAD::pow(vars(i, kIndexDelta) - vars(i - 1, kIndexDelta), 2);
      fg_vector[0] += kCostDiffA * CppAD::pow(vars(i, kIndexA) - vars(i - 1, kIndexA), 2);
    }

    fg_vector.segment(1, kExStateSize) = vars.row(0);

    for (int i = 0; i != fg.rows(); ++i) {
      auto x1 = vars(i + 1, kIndexX);
      auto y1 = vars(i + 1, kIndexY);
      auto psi1 = vars(i + 1, kIndexPsi);
      auto v1 = vars(i + 1, kIndexV);
      auto cte1 = vars(i + 1, kIndexCte);
      auto epsi1 = vars(i + 1, kIndexEpsi);

      auto x0 = vars(i, kIndexX);
      auto y0 = vars(i, kIndexY);
      auto psi0 = vars(i, kIndexPsi);
      auto v0 = vars(i, kIndexV);
      auto cte0 = vars(i, kIndexCte);
      auto epsi0 = vars(i, kIndexEpsi);

      auto delta0 = vars(i + 1, kIndexDelta);
      auto a0 = vars(i + 1, kIndexA);

      auto f0 = poly_(x0);
      auto psides0 = CppAD::atan(poly_(x0, 1));

      fg(i, kIndexX) = x1 - (x0 + v0 * CppAD::cos(psi0) * kTimeDelta);
      fg(i, kIndexY) = y1 - (y0 + v0 * CppAD::sin(psi0) * kTimeDelta);
      fg(i, kIndexPsi) = psi1 - (psi0 + v0 / kLf * delta0 * kTimeDelta);
      fg(i, kIndexV) = v1 - (v0 + a0 * kTimeDelta);
      fg(i, kIndexCte) = cte1 - (f0 - y0 + v0 * CppAD::sin(epsi0) * kTimeDelta);
      fg(i, kIndexEpsi)= epsi1 - (psi0 - psides0 + v0 / kLf * delta0 * kTimeDelta);
    }
  }

private:
  Poly<CppAD::AD<double>> poly_;
  double speed_;
};

Eigen::MatrixXd Mpc::Solve(double delay, double x0, double y0, double psi0, double v0, double delta, double a,
                           const vector<double>& ptsx,
                           const vector<double>& ptsy) {
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

  auto poly = Poly<double>::Fit(convert<VectorXd>(ptsx), convert<VectorXd>(ptsy), 3);

  double epsi0 = psi0 - atan(poly(x0, 1));

  double x = x0 + v0 * cos(psi0) * delay;
  double y = y0 + v0 * sin(psi0) * delay;
  double psi = psi0 + v0 / kLf * delta *delay;
  double v = v0 + a * delay;
  double cte = poly(x0) - y0 + v0 * sin(epsi0) * delay;
  double epsi = psi0 - atan(poly(x0, 1)) + v0 / kLf * delta * delay;

  VectorXd vars(kSteps * kExStateSize);
  vars.setZero();
  vars.segment(0, kExStateSize)
    << x, y, psi, v, cte, epsi, delta, a;

  VectorXd vars_lowerbound(vars.size());
  VectorXd vars_upperbound(vars.size());

  for (size_t i = 0; i < kSteps; ++i) {
    vars_lowerbound.segment(i * kExStateSize, kStateSize).array() = -1.0e19;
    vars_upperbound.segment(i * kExStateSize, kStateSize).array() = 1.0e19;

    vars_lowerbound[i * kExStateSize + kIndexDelta] = -0.436332;
    vars_upperbound[i * kExStateSize + kIndexDelta] = 0.436332;

    vars_lowerbound[i * kExStateSize + kIndexA] = -1;
    vars_upperbound[i * kExStateSize + kIndexA] = 1;
  }

  VectorXd constraints_lowerbound(kSteps * kStateSize + 2);
  VectorXd constraints_upperbound(kSteps * kStateSize + 2);
  constraints_lowerbound.setZero();
  constraints_upperbound.setZero();
  constraints_lowerbound.head(kExStateSize) = vars.head(kExStateSize);
  constraints_upperbound.head(kExStateSize) = vars.head(kExStateSize);

  FgEval fg_eval(poly, speed_ * speed_coeff_);

  std::string options;
  options += "Integer print_level  0\n";
  options += "String  sb           yes\n";
  options += "Sparse  true         forward\n";
  options += "Sparse  true         reverse\n";
  options += "Numeric max_cpu_time 0.5\n";

  CppAD::ipopt::solve_result<VectorXd> solution;
  CppAD::ipopt::solve(options, vars, vars_lowerbound, vars_upperbound,
                      constraints_lowerbound, constraints_upperbound, fg_eval, solution);

  if (solution.status != CppAD::ipopt::solve_result<VectorXd>::success) {
    std::cerr << "Couldn't find solution" << std::endl;
    return Eigen::MatrixXd();
  }

  std::cout << "Cost " << solution.obj_value << std::endl;
  
  Eigen::Map<const Matrix> result(solution.x.data(), kSteps, kExStateSize);

  speed_coeff_ = 1. - std::max(fabs(result.col(kIndexDelta).minCoeff()),
                               fabs(result.col(kIndexDelta).maxCoeff()));
  
  return result;
}
