#ifndef MPC_MPC_H_
#define MPC_MPC_H_

#include <vector>

#include "Eigen/Core"


class Mpc {
 public:
  Mpc(double speed)
    : speed_(speed) {}
  
  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  Eigen::MatrixXd Solve(double delay, double x, double y, double psi, double v, double delta, double a,
                        const std::vector<double>& ptsx,
                        const std::vector<double>& ptsy);

private:
  double speed_ = 40;
  double speed_coeff_ = 1;
};

#endif // MPC_MPC_H_
