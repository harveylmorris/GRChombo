/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        // morris:
        // added in lambda and eta
        double scalar_mass;
        double lambda;
        double eta;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // morris: OLD:
        // The potential value at phi
        // 1/2 m^2 phi^2
        //V_of_phi = 0.5 * pow(m_params.scalar_mass * vars.phi, 2.0);
        // The potential gradient at phi
        // m^2 phi
        // dVdphi = pow(m_params.scalar_mass, 2.0) * vars.phi;

        // morris NEW:
        // The potential value at phi
        // lambda / 4 * (phi^2 - eta^2)^2
        // add in for 1,2,3
        V_of_phi = 0.25 * m_params.lambda * pow(pow(vars.phi, 2) - pow(m_params.eta, 2), 2);
        // The potential gradient at phi
        // lambda * phi * (phi^2 - eta^2)
        dVdphi = m_params.lambda * vars.phi * (pow(vars.phi, 2) - pow(m_params.eta, 2));
    }
};

#endif /* POTENTIAL_HPP_ */
