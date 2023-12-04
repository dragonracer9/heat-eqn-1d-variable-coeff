#include "create_poisson_matrix.hpp"

//! Used for filling the sparse matrix.
using Triplet = Eigen::Triplet<double>;

//! Create the 1D Poisson matrix
//! @param[in] N the number of interior points
//! @param[in] a the coefficient function a
//!
//! @returns the Poisson matrix.
SparseMatrix createPoissonMatrix(int N,
                                 const std::function<double(double)> &a) {
  SparseMatrix A;
  // (write your solution here)
  A.resize(N, N);
  const double h = 1.0 / (N + 1);
  const auto xs = Eigen::VectorXd::LinSpaced(
      N + 2, 0, 1);  // linspace(low,high,size)' (from the cheatsheet)

  std::vector<Triplet> triplets;
  triplets.reserve((std::size_t)(N + N * 2 - 2));

  for (int i{0}; i < N; ++i) {
    // diagonal
    triplets.push_back(Triplet(i, i, 2.0 * a(xs(i + 1)) / (h * h)));
    // off-diagonal (upper)
    if (i > 0) {
      triplets.push_back(Triplet(i, i - 1, -a(xs(i + 1)) / (h * h)));
    }
    // off-diagonal (lower)
    if (i < N - 1) {
      triplets.push_back(Triplet(i, i + 1, -a(xs(i + 1)) / (h * h)));
    }
  }
  A.setFromTriplets(triplets.begin(), triplets.end());
  return A;
}
