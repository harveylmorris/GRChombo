/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"
#include "DiagnosticVariables.hpp"

// assign an enum to each variable
enum
{
    // Note that it is important that the first enum value is set to 1 more than
    // the last CCZ4 var enum

    /////////////////////////////////////
    // OLD
    //c_phi = NUM_CCZ4_VARS, // matter field added
    //c_Pi,                  //(minus) conjugate momentum
    /////////////////////////////////////

    // NEW
    c_phi1 = NUM_CCZ4_VARS, // matter field 1 added
    c_phi2 = NUM_CCZ4_VARS, // matter field 2 added
    c_phi3 = NUM_CCZ4_VARS, // matter field 3 added
    c_Pi1,                  //(minus) conjugate momentum 1
    c_Pi2,                  //(minus) conjugate momentum 2
    c_Pi3,                  //(minus) conjugate momentum 3

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    /////////////////////////////////////
    // OLD
    //user_variable_names = {"phi", "Pi"};
    /////////////////////////////////////

    // NEW
    user_variable_names = {"phi1", "phi2", "phi3", "Pi1", "Pi2", "Pi3"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
