/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "InitialScalarData.hpp"
#include "KerrBH.hpp"
#include "Potential.hpp"

// morris: reading initial data from text file and storing contents in an array
// files phi.txt
#include <fstream>
#include <iostream>
#include <vector>

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial scalar field data
        initial_params.center =
            center; // already read in SimulationParametersBase
        pp.load("G_Newton", G_Newton,
                0.0); // for now the example neglects backreaction
        pp.load("scalar_amplitude", initial_params.amplitude, 0.1);
        pp.load("scalar_width", initial_params.width, 1.0);
        pp.load("scalar_mass", potential_params.scalar_mass, 0.1);

        // Initial Kerr data
        pp.load("kerr_mass", kerr_params.mass, 1.0);
        pp.load("kerr_spin", kerr_params.spin, 0.0);
        pp.load("kerr_center", kerr_params.center, center);

        // morris: adding parameters for potential
        pp.load("pot_lambda", potential_params.pot_lambda);
        pp.load("pot_eta", potential_params.pot_eta);
        initial_params.pot_eta = potential_params.pot_eta;

        int num_elements;
        pp.load("num_elements", num_elements);
        pp.load("spacing", initial_params.spacing);

	double twist_0;
        pp.load("twist", twist_0);
	initial_params.twist = 2.0 * M_PI * twist_0;
        pp.load("vel_z", initial_params.vel_z);

        pp.load("threshold_phi", threshold_phi);
        pp.load("threshold_K", threshold_K);

        double initial_f[num_elements];
        double initial_A[num_elements];

        ifstream initial_f_file("gr_initial_f_eta7e-2.txt");
        double tmp1 = 0;
        for (int i = 0; i < num_elements; ++i)
        {
            initial_f_file >> tmp1;
            initial_f[i] = tmp1;
        }
        initial_f_file.close();

        ifstream initial_A_file("gr_initial_A_eta7e-2.txt");
        double tmp2 = 1;
        for (int i = 0; i < num_elements; ++i)
        {
            initial_A_file >> tmp2;
            initial_A[i] = tmp2;
        }
        initial_A_file.close();

        initial_params.p_initial_f = initial_f;
        initial_params.p_initial_A = initial_A;

    }

    void check_params()
    {
        warn_parameter("scalar_mass", potential_params.scalar_mass,
                       potential_params.scalar_mass <
                           0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of scalar field do not appear to be "
                       "resolved on coarsest level");
        warn_parameter("scalar_width", initial_params.width,
                       initial_params.width < 0.5 * L,
                       "is greater than half the domain size");
        warn_parameter("kerr_mass", kerr_params.mass, kerr_params.mass >= 0.0,
                       "should be >= 0.0");
        check_parameter("kerr_spin", kerr_params.spin,
                        std::abs(kerr_params.spin) <= kerr_params.mass,
                        "must satisfy |a| <= M = " +
                            std::to_string(kerr_params.mass));
        FOR(idir)
        {
            std::string name = "kerr_center[" + std::to_string(idir) + "]";
            warn_parameter(
                name, kerr_params.center[idir],
                (kerr_params.center[idir] >= 0) &&
                    (kerr_params.center[idir] <= (ivN[idir] + 1) * coarsest_dx),
                "should be within the computational domain");
        }
    }

    // Initial data for matter and potential and BH
    double G_Newton, threshold_phi, threshold_K;
    InitialScalarData::params_t initial_params;
    Potential::params_t potential_params;
    KerrBH::params_t kerr_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
