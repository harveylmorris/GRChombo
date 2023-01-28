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
        data_t rr2 = rr * rr;

        // calculate the field value
        data_t phi = m_params.amplitude *
                     (1.0 + 0.01 * rr2 * exp(-pow(rr / m_params.width, 2.0)));

        // store the vars
        // TODO(morris): remove and add c_phi1, c_phi2, c_phi3, c_pi1, c_pi2, c_pi3
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
