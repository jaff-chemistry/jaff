// ABOUTME: Network header of the JAFF-generated chemistry for AMReX Microphysics:
// ABOUTME: species come from pynucastro.net, parameters from _parameters (network_rp)
#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <string_view>

#include <fundamental_constants.H>
#include <network_properties.H>
#include <burn_type.H>

// Checks that every network parameter is set (stops naming those that are not).
void actual_network_init();

constexpr std::string_view network_name = "jaff-chem";

namespace network
{
}

// Sets the electron number density from the charges of the other species.
void balance_charge(burn_t& state);

#endif
