"""
    HeatPipeCore

A one-dimensional, non-CFD heat-pipe simulator for space fission power
systems.  The package couples:

* alkali-metal thermophysical properties (`FluidProperties`);
* closed-form wick geometry and permeability (`Wick`);
* analytical operational limits (`Limits`);
* a lumped-parameter axial thermal resistance network (`ThermalNetwork`);
* benchmark definitions (`Validation`);
* the ISR Mark 2 design-space study (`ISRDesign`).

Typical use:

```julia
using HeatPipeCore
hp = krusty_heatpipe()
res = steady_state(hp, 300.0)
res.T_v - 273.15          # vapour temperature in degrees Celsius
res.limits.minimum        # controlling power limit at that temperature
max_power(hp, 1073.15)    # power limit at 800 degC
```
"""
module HeatPipeCore

include("FluidProperties.jl")
include("Wick.jl")
include("Limits.jl")
include("ThermalNetwork.jl")
include("Validation.jl")
include("ISRDesign.jl")

using .FluidProperties
using .Wick
using .Limits
using .ThermalNetwork
using .Validation
using .ISRDesign

export Fluid, Na, K, Li, fluid, p_sat, rho_l, rho_v, h_fg, sigma,
       mu_l, mu_v, k_l, cp_l, T_melt, T_boil, T_crit, R_specific,
       WickGeometry, screen_wick, sintered_wick, annular_wick,
       k_eff_maxwell, HeatPipe, q_viscous, q_sonic, q_entrainment,
       q_capillary, q_boiling, q_flooding, operating_limits, max_power,
       controlling_limit, effective_length, vapour_reynolds,
       ThermalNetworkModel, build_network, steady_state, transient_response,
       stable_dt, krusty_heatpipe, krusty_benchmark, krusty_validation,
       krusty_limit, krusty_flood_empirical, walker_heatpipe,
       walker_validation, nominal_rows, nominal_point_validation,
       validation_table,
       isr_core, isr_heatpipe, isr_pin_power, isr_feasibility,
       isr_diameter_sweep, isr_fluid_comparison, isr_options,
       isr_mass_balance, isr_required_diameter

end # module
