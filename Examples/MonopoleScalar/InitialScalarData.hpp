/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"


// morris: code for linear interpolation
#include <iostream>
#include <vector>

double linearInterpolation(double x, const std::vector<double> &X, const std::vector<double> &Y) {
  int n = X.size();
  if (n != Y.size()) {
    // X and Y should have the same size
    std::cerr << "Error: X and Y must have the same size." << std::endl;
    return 0;
  }

  if (x <= X[0]) {
    // If x is less than or equal to the first value of X, return the first value of Y
    return Y[0];
  }

  if (x >= X[n - 1]) {
    // If x is greater than or equal to the last value of X, return the last value of Y
    return Y[n - 1];
  }

  // Find the index of the largest value in X that is less than or equal to x
  auto it = std::lower_bound(X.begin(), X.end(), x);
  int i = it - X.begin() - 1;

  // Perform linear interpolation between Y[i] and Y[i + 1]
  return Y[i] + (Y[i + 1] - Y[i]) * (x - X[i]) / (X[i + 1] - X[i]);
}


//! Class which sets the initial scalar field matter config
class InitialScalarData
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //!< Amplitude of bump in initial SF bubble
        std::array<double, CH_SPACEDIM>
            center;   //!< Centre of perturbation in initial SF bubble
        double width; //!< Width of bump in initial SF bubble
        // morris
        double eta;
        vector<float> initial_f;
        vector<float> initial_r;
    };

    //! The constructor
    InitialScalarData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        data_t rr = coords.get_radius();

        // morris: calculate the field value
        // field configuration describing a monopole is phi^a = eta * f(r) * x^a / r
        // first we find f based on r
        double f = linearInterpolation(rr, m_params.initial_r, m_params.initial_f)
        // then we find phi^a
        // starting with just one scalar field so let's use x coordinate
        // JOSU: should I used coords.x or m_params.center[0]?
        if (rr == 0) {
            data_t phi = m_params.eta * f;
        } else {
            data_t phi = m_params.eta * f * coords.x / rr;
        }
        
        // store the vars
        // TODO(morris): remove // and switch to 3 dimensions
        current_cell.store_vars(phi, c_phi);
        current_cell.store_vars(0.0, c_Pi);
        // current_cell.store_vars(phi1, c_phi1);
        // current_cell.store_vars(phi2, c_phi2);
        // current_cell.store_vars(phi3, c_phi3);
        // current_cell.store_vars(0.0, c_Pi1);
        // current_cell.store_vars(0.0, c_Pi2);
        // current_cell.store_vars(0.0, c_Pi3);

        // morris: adding metric components
        current_cell.store_vars(1.0, c_h11);
        current_cell.store_vars(1.0, c_h22);
        current_cell.store_vars(1.0, c_h33);
        
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALSCALARDATA_HPP_ */





