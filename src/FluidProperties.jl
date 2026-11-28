"""
    FluidProperties

Thermophysical properties of the alkali-metal working fluids used in the
simulator: sodium (Na), potassium (K) and lithium (Li).

All quantities are expressed in pure SI units:

| quantity                    | symbol   | unit    |
|:----------------------------|:---------|:--------|
| temperature                 | `T`      | K       |
| pressure                    | `p`      | Pa      |
| density                     | `rho`    | kg/m^3  |
| latent heat of vaporisation | `h_fg`   | J/kg    |
| surface tension             | `sigma`  | N/m     |
| dynamic viscosity           | `mu`     | Pa s    |
| thermal conductivity        | `k`      | W/(m K) |
| specific heat capacity      | `cp`     | J/(kg K)|

The saturation pressure is represented by a two-constant Clausius-Clapeyron
fit, `log10(p) = A - B/T`, fitted to tabulated alkali-metal vapour-pressure
data over the range of interest (roughly 600-1600 K).  The remaining
properties use standard engineering correlations:

* liquid density: linear in temperature;
* latent heat: Watson correlation, `h_fg = h_fg,b ((Tc-T)/(Tc-Tb))^0.38`;
* surface tension: Watson-type power law;
* liquid viscosity: Andrade (Arrhenius) form, `mu = mu0 exp(mu_B/T)`;
* vapour density: ideal gas, `rho_v = p M / (R T)`;
* vapour viscosity: weak power law in temperature.

Vapour-phase non-ideality is neglected, which is acceptable at the
moderate reduced pressures (p/p_c < 0.1) considered here.
"""
module FluidProperties

export Fluid, Na, K, Li, fluid, p_sat, rho_l, rho_v, h_fg, sigma,
       mu_l, mu_v, k_l, cp_l, T_melt, T_boil, T_crit, molar_mass,
       gamma_v, R_specific, dprops, is_valid

const R_UNIVERSAL = 8.31446261815324  # J/(mol K)

"""
    Fluid

Immutable description of an alkali-metal working fluid.
"""
struct Fluid
    name::String
    Mw::Float64        # molar mass, kg/mol
    T_melt::Float64    # melting point at 1 atm, K
    T_boil::Float64    # normal boiling point, K
    T_crit::Float64    # critical temperature, K
    logP_A::Float64    # log10(p/Pa) = A - B/T
    logP_B::Float64    # K
    h_fg_boil::Float64 # latent heat at the normal boiling point, J/kg
    rho0::Float64      # liquid density at T_melt, kg/m^3
    rho_slope::Float64 # d(rho_l)/dT, kg/(m^3 K), positive
    sigma0::Float64    # surface tension at T_melt, N/m
    mu0::Float64       # liquid viscosity prefactor, Pa s
    mu_B::Float64      # liquid viscosity Arrhenius constant, K
    mu_v0::Float64     # vapour viscosity at 1000 K, Pa s
    k_l::Float64       # liquid thermal conductivity, W/(m K)
    cp_l::Float64      # liquid specific heat, J/(kg K)
    gamma::Float64     # ratio of specific heats of the vapour
end

"""
    Na, K, Li

Pre-defined working fluids.
"""
const Na = Fluid("Sodium", 0.0229898, 370.94, 1156.1, 2573.0,
                 9.710, 5438.0, 4.24e6,
                 927.0, 0.235, 0.206,
                 1.20e-4, 700.0, 1.50e-5,
                 72.0, 1270.0, 5.0 / 3.0)

const K = Fluid("Potassium", 0.0390983, 336.7, 1032.0, 2223.0,
                9.408, 4543.0, 1.95e6,
                828.0, 0.270, 0.110,
                1.00e-4, 800.0, 1.30e-5,
                45.0, 760.0, 5.0 / 3.0)

const Li = Fluid("Lithium", 0.006941, 453.7, 1615.0, 3223.0,
                 10.061, 8164.0, 20.9e6,
                 535.0, 0.105, 0.400,
                 3.00e-5, 1200.0, 1.20e-5,
                 45.0, 4200.0, 5.0 / 3.0)

const _FLUIDS = Dict("Na" => Na, "K" => K, "Li" => Li, "Sodium" => Na,
                     "Potassium" => K, "Lithium" => Li)

"""Return a fluid by symbol or name, e.g. `fluid(:Na)`."""
fluid(x::Fluid) = x
fluid(x::Symbol) = _FLUIDS[String(x)]
fluid(x::AbstractString) = _FLUIDS[x]

molar_mass(f::Fluid) = f.Mw
R_specific(f::Fluid) = R_UNIVERSAL / f.Mw
gamma_v(f::Fluid) = f.gamma
T_melt(f::Fluid) = f.T_melt
T_boil(f::Fluid) = f.T_boil
T_crit(f::Fluid) = f.T_crit

"""Saturation pressure `p_sat(T)` in Pa (`log10 p = A - B/T`)."""
p_sat(f::Fluid, T::Real) = 10.0^(f.logP_A - f.logP_B / T)

"""Saturated liquid density in kg/m^3 (linear correlation)."""
rho_l(f::Fluid, T::Real) = max(f.rho0 - f.rho_slope * (T - f.T_melt), 0.05 * f.rho0)

"""Saturated vapour density in kg/m^3 from the ideal-gas law."""
rho_v(f::Fluid, T::Real) = p_sat(f, T) * f.Mw / (R_UNIVERSAL * T)

"""Latent heat of vaporisation in J/kg (Watson correlation)."""
function h_fg(f::Fluid, T::Real)
    T = clamp(float(T), f.T_melt, 0.999 * f.T_crit)
    return f.h_fg_boil * ((f.T_crit - T) / (f.T_crit - f.T_boil))^0.38
end

"""Surface tension in N/m (Watson-type power law)."""
function sigma(f::Fluid, T::Real)
    T = clamp(float(T), f.T_melt, 0.999 * f.T_crit)
    return f.sigma0 * ((f.T_crit - T) / (f.T_crit - f.T_melt))^1.22
end

"""Saturated liquid dynamic viscosity in Pa s (Andrade form)."""
mu_l(f::Fluid, T::Real) = f.mu0 * exp(f.mu_B / T)

"""Saturated vapour dynamic viscosity in Pa s."""
mu_v(f::Fluid, T::Real) = f.mu_v0 * (T / 1000.0)^0.7

"""Liquid thermal conductivity in W/(m K)."""
k_l(f::Fluid, T::Real) = f.k_l

"""Liquid specific heat capacity in J/(kg K)."""
cp_l(f::Fluid, T::Real) = f.cp_l

"""Return `true` when `T` lies inside the liquid-vapour coexistence range."""
is_valid(f::Fluid, T::Real) = f.T_melt <= T < f.T_crit

"""
    dprops(f, T) -> NamedTuple

Bundle of all temperature-dependent properties at saturation, convenient for
logging and for the report tables.
"""
function dprops(f::Fluid, T::Real)
    return (T = float(T),
            p_sat = p_sat(f, T),
            rho_l = rho_l(f, T),
            rho_v = rho_v(f, T),
            h_fg = h_fg(f, T),
            sigma = sigma(f, T),
            mu_l = mu_l(f, T),
            mu_v = mu_v(f, T),
            k_l = k_l(f, T),
            cp_l = cp_l(f, T))
end

end # module
