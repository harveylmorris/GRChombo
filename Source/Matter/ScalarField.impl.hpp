/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SCALARFIELD_HPP_)
#error "This file should only be included through ScalarField.hpp"
#endif

#ifndef SCALARFIELD_IMPL_HPP_
#define SCALARFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> ScalarField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    // call the function which computes the em tensor excluding the potential
    emtensor_excl_potential(out, vars, d1, h_UU, chris_ULL);

    /////////////////////////////////////
    // OLD
    // set the potential values
    //data_t V_of_phi = 0.0;
    //data_t dVdphi = 0.0;

    // compute potential and add constributions to EM Tensor
    //my_potential.compute_potential(V_of_phi, dVdphi, vars);

    //out.rho += V_of_phi;
    //out.S += -3.0 * V_of_phi;

    // FOR(i, j) { out.Sij[i][j] += -vars.h[i][j] * V_of_phi / vars.chi; }
    /////////////////////////////////////

    // morris: NEW with 3 scalar fields:
    // set the potential values
    data_t V_of_phi1 = 0.0;
    data_t V_of_phi2 = 0.0;
    data_t V_of_phi3 = 0.0;
    data_t dVdphi1 = 0.0;
    data_t dVdphi2 = 0.0;
    data_t dVdphi3 = 0.0;

    // compute potential and add constributions to EM Tensor
    my_potential.compute_potential(V_of_phi1, V_of_phi2, V_of_phi3,
                                   dVdphi1, dVdphi2, dVdphi3, vars);

    out.rho += V_of_phi1 + V_of_phi2 + V_of_phi3;
    out.S += -3.0 * (V_of_phi1 + V_of_phi2 + V_of_phi3);

    FOR(i, j) { out.Sij[i][j] += -vars.h[i][j] * (V_of_phi1 + V_of_phi2 + V_of_phi3) / vars.chi; }

    return out;
}

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void ScalarField<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &h_UU,
    const Tensor<3, data_t> &chris_ULL)
{
    /////////////////////////////////////
    // OLD
    // Useful quantity Vt
    //data_t Vt = -vars.Pi * vars.Pi;
    //FOR(i, j) { Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    //FOR(i, j)
    //{
    //    out.Sij[i][j] =
    //        -0.5 * vars.h[i][j] * Vt / vars.chi + d1.phi[i] * d1.phi[j];
    //}

    // S = Tr_S_ij
    //out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // S_i (note lower index) = - n^a T_ai
    //FOR(i) { out.Si[i] = -d1.phi[i] * vars.Pi; }

    // rho = n^a n^b T_ab
    //out.rho = vars.Pi * vars.Pi + 0.5 * Vt;
    /////////////////////////////////////


    // NEW
    // Useful quantity Vt
    data_t Vt = -vars.Pi1 * vars.Pi1 -vars.Pi2 * vars.Pi2 -vars.Pi3 * vars.Pi3;
    FOR(i, j) { Vt += vars.chi * h_UU[i][j] * (d1.phi1[i] * d1.phi1[j] + d1.phi2[i] * d1.phi2[j] + d1.phi3[i] * d1.phi3[j]); }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR(i, j)
    {
        out.Sij[i][j] =
            -0.5 * vars.h[i][j] * Vt / vars.chi + (d1.phi1[i] * d1.phi1[j] + d1.phi2[i] * d1.phi2[j] + d1.phi3[i] * d1.phi3[j]);
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // S_i (note lower index) = - n^a T_ai
    FOR(i) { out.Si[i] = -(d1.phi1[i] * vars.Pi1 + d1.phi2[i] * vars.Pi2 + d1.phi3[i] * vars.Pi3); }

    // rho = n^a n^b T_ab
    out.rho = (vars.Pi1 * vars.Pi1 + vars.Pi2 * vars.Pi2 + vars.Pi3 * vars.Pi3) + 0.5 * Vt;

}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void ScalarField<potential_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // first get the non potential part of the rhs
    // this may seem a bit long winded, but it makes the function
    // work for more multiple fields

    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(total_rhs, vars, d1, d2, advec);

    /////////////////////////////////////
    // OLD
    // set the potential values
    //data_t V_of_phi = 0.0;
    //data_t dVdphi = 0.0;
    //my_potential.compute_potential(V_of_phi, dVdphi, vars);

    // adjust RHS for the potential term
    //total_rhs.Pi += -vars.lapse * dVdphi;
    /////////////////////////////////////

    // NEW
    data_t V_of_phi1 = 0.0;
    data_t V_of_phi2 = 0.0;
    data_t V_of_phi3 = 0.0;
    data_t dVdphi1 = 0.0;
    data_t dVdphi2 = 0.0;
    data_t dVdphi3 = 0.0;
    my_potential.compute_potential(V_of_phi1, V_of_phi2, V_of_phi3,
                                   dVdphi1, dVdphi2, dVdphi3, vars);
    // adjust RHS for the potential term
    total_rhs.Pi1 += -vars.lapse * dVdphi1;
    total_rhs.Pi2 += -vars.lapse * dVdphi2;
    total_rhs.Pi3 += -vars.lapse * dVdphi3;

}

// the RHS excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void ScalarField<potential_t>::matter_rhs_excl_potential(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec)
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // evolution equations for scalar field and (minus) its conjugate momentum

    /////////////////////////////////////
    // OLD
    //rhs.phi = vars.lapse * vars.Pi + advec.phi;
    //rhs.Pi = vars.lapse * vars.K * vars.Pi + advec.Pi;

    //FOR(i, j)
    //{ 
        // includes non conformal parts of chris not included in chris_ULL
    //    rhs.Pi += h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1.phi[i] +
    //                            vars.chi * vars.lapse * d2.phi[i][j] +
    //                            vars.chi * d1.lapse[i] * d1.phi[j]);
    //    FOR(k)
    //    {
    //        rhs.Pi += -vars.chi * vars.lapse * h_UU[i][j] * chris.ULL[k][i][j] *
    //                  d1.phi[k];
    //    }
    //}
    /////////////////////////////////////
    

    // NEW
    rhs.phi1 = vars.lapse * vars.Pi1 + advec.phi1;
    rhs.phi2 = vars.lapse * vars.Pi2 + advec.phi2;
    rhs.phi3 = vars.lapse * vars.Pi3 + advec.phi3;
    rhs.Pi1 = vars.lapse * vars.K * vars.Pi1 + advec.Pi1;
    rhs.Pi2 = vars.lapse * vars.K * vars.Pi2 + advec.Pi2;
    rhs.Pi3 = vars.lapse * vars.K * vars.Pi3 + advec.Pi3;

    FOR(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs.Pi1 += h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1.phi1[i] +
                                vars.chi * vars.lapse * d2.phi1[i][j] +
                                vars.chi * d1.lapse[i] * d1.phi1[j]);
        rhs.Pi2 += h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1.phi2[i] +
                                vars.chi * vars.lapse * d2.phi2[i][j] +
                                vars.chi * d1.lapse[i] * d1.phi2[j]);
        rhs.Pi3 += h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1.phi3[i] +
                                vars.chi * vars.lapse * d2.phi3[i][j] +
                                vars.chi * d1.lapse[i] * d1.phi3[j]);
        FOR(k)
        {
            rhs.Pi1 += -vars.chi * vars.lapse * h_UU[i][j] * chris.ULL[k][i][j] *
                      d1.phi1[k];
            rhs.Pi2 += -vars.chi * vars.lapse * h_UU[i][j] * chris.ULL[k][i][j] *
                      d1.phi2[k];
            rhs.Pi3 += -vars.chi * vars.lapse * h_UU[i][j] * chris.ULL[k][i][j] *
                      d1.phi3[k];
        }
    }
}

#endif /* SCALARFIELD_IMPL_HPP_ */
