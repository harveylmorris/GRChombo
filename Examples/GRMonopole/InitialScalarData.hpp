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
        double *p_initial_A;
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

        double rho = sqrt(coords.x * coords.x + coords.y * coords.y +
                           coords.z * coords.z);
        
        double xx = coords.x;
        double yy = coords.y;
        double zz = coords.z;

        data_t rr = coords.get_radius();

        // field configuration describing a monopole is phi^a = eta * f(r) * x^a
        // / r first we find f based on r

        int indxL = static_cast<int>(floor(rho / m_params.spacing));
        int indxH = static_cast<int>(ceil(rho / m_params.spacing));
        double f_data_L = *(m_params.p_initial_f + indxL);
        double f_data_H = *(m_params.p_initial_f + indxH);
        double A_data_L = *(m_params.p_initial_A + indxL);
        double A_data_H = *(m_params.p_initial_A + indxH);

        data_t f = f_data_L + (rho / m_params.spacing - indxL) * (f_data_H - f_data_L);
        data_t A = A_data_L + (rho / m_params.spacing - indxL) * (A_data_H - A_data_L);
        data_t B = 1 / A;

        data_t phi1 = m_params.pot_eta * f * coords.x / rr;
        data_t phi2 = m_params.pot_eta * f * coords.y / rr;
        data_t phi3 = m_params.pot_eta * f * coords.z / rr;

        current_cell.store_vars(phi1, c_phi1);
        current_cell.store_vars(phi2, c_phi2);
        current_cell.store_vars(phi3, c_phi3);

        current_cell.store_vars(0.0, c_Pi1);
        current_cell.store_vars(0.0, c_Pi2);
        current_cell.store_vars(0.0, c_Pi3);

        // Adding metric components for GR:

        double gamma_sph[3][3];

        // metric: ds^2 = B(r) dt^2 - A(r) dr^2 - r^2 (d{theta}^2 + sin^2(theta) d{phi}^2)
        // so:
        // gamma_rr = A(r)
        // gamma_thetatheta = r^2
        // gamma_phiphi = r^2 sin^2(theta)
        // sin^2(theta) = 1 - cos^2(theta) = 1 - (z/r)^2

        gamma_sph[0][0] = A;
        gamma_sph[1][1] = rho * rho;
        gamma_sph[2][2] = rho * rho * (1 - pow(zz / rho, 2));

        // Define jacobian for change of coordinates

        double jacobian[3][3];

        jacobian[0][0] = xx/rho;
        jacobian[1][0] = xx * zz / (pow(rho, 3) * pow(1 - pow(zz/rho, 2), 0.5));
        jacobian[2][0] = - yy / (1 - pow(zz/rho, 2));
        
        jacobian[0][1] = yy/rho;
        jacobian[1][1] = yy * zz / (pow(rho, 3) * pow(1 - pow(zz/rho, 2), 0.5));
        jacobian[2][1] = xx / (1 - pow(zz/rho, 2));
        
        jacobian[0][2] = zz/rho;
        jacobian[1][2] = - pow(1 - pow(zz/rho, 2), 0.5) / rho;
        jacobian[2][2] = 0.0;

        // Coordinate transformation : Spherical -> Cartesian
        double gammaxyz[3][3];

        FOR2(i,j)
        {
            FOR2(k,l){
                gammaxyz[i][j] += gamma_sph[k][l]*jacobian[k][i]*jacobian[l][j];
            }
        }

        // conformal decomposition

        data_t deth = gammaxyz[0][0] * (gammaxyz[1][1]*gammaxyz[2][2] - gammaxyz[1][2]*gammaxyz[2][1]) -
                    gammaxyz[0][1] * (gammaxyz[2][2]*gammaxyz[1][0] - gammaxyz[1][2]*gammaxyz[2][0]) + 
                    gammaxyz[0][2] * (gammaxyz[1][0]*gammaxyz[2][1] - gammaxyz[1][1]*gammaxyz[2][0]);

        data_t chi = pow(deth, -1./3.);

        // QUESTION: not sure on lapse?
        current_cell.store_vars(1.0, c_lapse); 
        current_cell.store_vars(chi, c_chi);

        current_cell.store_vars(gammaxyz[0][0], c_h11);
        current_cell.store_vars(gammaxyz[0][1], c_h12);
        current_cell.store_vars(gammaxyz[0][2], c_h13);

        current_cell.store_vars(gammaxyz[1][1], c_h22);
        current_cell.store_vars(gammaxyz[1][2], c_h23);
        
        current_cell.store_vars(gammaxyz[2][2], c_h33);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALSCALARDATA_HPP_ */
