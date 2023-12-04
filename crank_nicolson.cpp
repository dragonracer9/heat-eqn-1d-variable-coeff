#include "crank_nicolson.hpp"

//! Uses the Crank-Nicolson method to compute u from time 0 to time T
//!
//! @param[in] u0 the initial data, as column vector (size N+2)
//! @param[in] dt the time step size
//! @param[in] T the final time at which to compute the solution (which we assume to be a multiple of dt)
//! @param[in] N the number of interior grid points
//! @param[in] gL function of time with the Dirichlet condition at left boundary
//! @param[in] gR function of time with the Dirichlet condition at right boundary
//! @param[in] a the coefficient function a
//!
//! @note the vector returned should include the boundary values!
//!
std::pair<Eigen::MatrixXd, Eigen::VectorXd> crankNicolson(
    const Eigen::VectorXd& u0,
    double dt, double T, int N,
    const std::function<double(double)>& gL,
    const std::function<double(double)>& gR,
    const std::function<double(double)>& a)
{

    Eigen::VectorXd time;
    Eigen::MatrixXd u;

    // (write your solution here)
    SparseMatrix I(N, N);
    I.setIdentity(); // shame this can't be done in the constructor (non-const)

    // create the matrix A
    SparseMatrix A = createPoissonMatrix(N, a);

    // set up u and time
    u.resize(N + 2, int(round(T / dt)) + 1);
    time.resize(int(round(T / dt)) + 1);

    // set up the initial conditions
    u.col(0) = u0;
    time(0) = 0.0;

    const auto h = 1. / (N + 1);

    // initialise the right hand side
    Eigen::VectorXd G_k_plus_one = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd G_k = Eigen::VectorXd::Zero(N);

    // set up solver
    Eigen::SparseLU<SparseMatrix> solver; // SimplicialLDLT<SparseMatrix> fails the last test for some reason
    solver.compute(I + 0.5 * dt * A);

    SparseMatrix B = I - 0.5 * dt * A;

    // iterate over time
    for (int i { 0 }; i < int(round(T / dt)); ++i) {

        // update time
        time(i + 1) = (i + 1) * dt;

        // set up the right hand side
        G_k(0) = dt * gL(time(i)) * a(h) / (h * h); // a(h) = a(x_1)
        G_k(N - 1) = dt * gR(time(i)) * a(1 - h) / (h * h); // a(1 - h) = a(x_N)

        G_k_plus_one(0) = dt * gL(time(i + 1)) * a(h) / (h * h); // a(h) = a(x_1)
        G_k_plus_one(N - 1) = dt * gR(time(i + 1)) * a(1 - h) / (h * h); // a(1 - h) = a(x_N)

        // solve the system
        u.col(i + 1).segment(1, N) = solver.solve(B * u.col(i).segment(1, N) + 0.5 * (G_k_plus_one + G_k)); // leave the endpoints alone

        // update the boundary conditions
        u.col(i + 1)(0) = gL(time(i + 1));
        u.col(i + 1)(N + 1) = gR(time(i + 1));
    }

    return std::make_pair(u, time);
}
