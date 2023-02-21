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

// for floor and ceil
#include <cmath>

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
        std::array<double, CH_SPACEDIM>
            center_monopole1;   //!< Centre of perturbation in initial SF bubble
        std::array<double, CH_SPACEDIM>
            center_monopole2;   //!< Centre of perturbation in initial SF bubble
        double width; //!< Width of bump in initial SF bubble
        // morris
        double pot_eta;
        double *p_initial_f;
        double *p_initial_f_prime;
        double spacing;
        double twist;
	double vel_z;
    };

    //! The constructor
    InitialScalarData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i relative to center?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        // where am i relative to the 2 monopoles?
        //Coordinates<data_t> coords1(current_cell, m_dx, m_params.center_monopole1);
        //data_t rr1 = coords1.get_radius();
        //Coordinates<data_t> coords2(current_cell, m_dx, m_params.center_monopole2);
        //data_t rr2 = coords2.get_radius();

        //double rho1 = sqrt(coords1.x * coords1.x + coords1.y * coords1.y +
        //                   coords1.z * coords1.z);
        //double rho2 = sqrt(coords2.x * coords2.x + coords2.y * coords2.y +
        //                   coords2.z * coords2.z);

        // field configuration describing a monopole is phi^a = eta * f(r) * x^a
        // / r first we find f based on r

        //int indxL1 = static_cast<int>(floor(rho1 / m_params.spacing));
        //int indxH1 = static_cast<int>(ceil(rho1 / m_params.spacing));
        //double f_data_L1 = *(m_params.p_initial_f + indxL1);
        //double f_data_H1 = *(m_params.p_initial_f + indxH1);

        //int indxL2 = static_cast<int>(floor(rho2 / m_params.spacing));
        //int indxH2 = static_cast<int>(ceil(rho2 / m_params.spacing));
        //double f_data_L2 = *(m_params.p_initial_f + indxL2);
        //double f_data_H2 = *(m_params.p_initial_f + indxH2);

        //data_t f1 =
        //    f_data_L1 + (rho1 / m_params.spacing - indxL1) * (f_data_H1 - f_data_L1);

        //data_t f2 =
        //    f_data_L2 + (rho2 / m_params.spacing - indxL2) * (f_data_H2 - f_data_L2);
        
        /////////////////////////////////////
        // OLD
        //data_t phi = m_params.pot_eta * f * coords.x / rr;
        // store the vars
        //current_cell.store_vars(phi, c_phi);
        //current_cell.store_vars(0.0, c_Pi);
        /////////////////////////////////////

        // NEW
        // eqns 22, 23, 24 from https://arxiv.org/pdf/1705.03091.pdf

        //double s = sin(m_params.twist);
        //double c = cos(m_params.twist);
        //double z_0 = sqrt((coords1.z - coords.z) * (coords1.z - coords.z));

        //data_t phi1 = m_params.pot_eta * f1 * f2 * ((c * coords.x + s * coords.y) *
        //                                   ((coords.z + z_0) * c - (coords.z - z_0))
        //                                   - (c * coords.y - s * coords.x) * rr2 * s) / rr1 / rr2;
        //data_t phi2 = m_params.pot_eta * f1 * f2 * ((c * coords.y - s * coords.x) *
        //                                   ((coords.z + z_0) * c - (coords.z - z_0))
        //                                   + (c * coords.x + s * coords.y) * rr2 * s) / rr1 / rr2;
        //data_t phi3 = m_params.pot_eta * f1 * f2 * ((coords.z - z_0) * (coords.z + z_0) +
        //                                  (coords.x * coords.x + coords.y * coords.y) * c) / rr1 / rr2;
        // store the vars

        // BACK TO ONE MONOPOLE CASE BUT NOW WITH GALIEAN BOOST

        data_t rr = coords.get_radius();

        double rho = sqrt(coords.x * coords.x + coords.y * coords.y +
                          coords.z * coords.z);

        // field configuration describing a monopole is phi^a = eta * f(r) * x^a
        // / r first we find f based on r

        int indxL = static_cast<int>(floor(rho / m_params.spacing));
        int indxH = static_cast<int>(ceil(rho / m_params.spacing));
        double f_data_L = *(m_params.p_initial_f + indxL);
        double f_data_H = *(m_params.p_initial_f + indxH);
        double f_prime_data_L = *(m_params.p_initial_f_prime + indxL);
        double f_prime_data_H = *(m_params.p_initial_f_prime + indxH);

        data_t f = f_data_L + (rho / m_params.spacing - indxL) * (f_data_H - f_data_L);
        data_t f_prime = f_prime_data_L + (rho / m_params.spacing - indxL) * (f_prime_data_H - f_prime_data_L);

        data_t phi1 = m_params.pot_eta * f * coords.x / rr;
        data_t phi2 = m_params.pot_eta * f * coords.y / rr;
        data_t phi3 = m_params.pot_eta * f * coords.z / rr;

        current_cell.store_vars(phi1, c_phi1);
        current_cell.store_vars(phi2, c_phi2);
        current_cell.store_vars(phi3, c_phi3);

        // boosting monopole
        //float v_x = 0;
        //float v_y = 0;
        //float v_z = 0;
        //float pi1 = v_x * d1.phi1[1];
        //float pi2 = v_y * d1.phi2[2];
        //float pi3 = v_z * d1.phi3[3];

        double v = m_params.vel_z;  // speed at which monopole is moving

        data_t pi1 = m_params.pot_eta * (- f_prime * (v * coords.z / rr) * (coords.x / rr) + f * (v * coords.z / rr) * (coords.x / rr / rr));
        data_t pi2 = m_params.pot_eta * (- f_prime * (v * coords.z / rr) * (coords.y / rr) + f * (v * coords.z / rr) * (coords.y / rr / rr));
        data_t pi3 = m_params.pot_eta * (- f_prime * (v * coords.z / rr) * (coords.z / rr) + f * (v * coords.z / rr) * (coords.z / rr / rr) - f * v / rr);

        current_cell.store_vars(pi1, c_Pi1);
        current_cell.store_vars(pi2, c_Pi2);
        current_cell.store_vars(pi3, c_Pi3);

        // morris: adding metric components
        current_cell.store_vars(1.0, c_lapse); 
        current_cell.store_vars(1.0, c_chi);
        current_cell.store_vars(1.0, c_h11);
        current_cell.store_vars(1.0, c_h22);
        current_cell.store_vars(1.0, c_h33);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALSCALARDATA_HPP_ */
