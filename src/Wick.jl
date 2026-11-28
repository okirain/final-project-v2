"""
    Wick

Closed-form geometry and transport correlations for the three wick families
considered in this study:

* `screen_wick`   - woven wire screen (square mesh);
* `sintered_wick` - sintered spherical powder;
* `annular_wick`  - annular gap (or open groove) with a fine-pore surface mesh.

For each wick the module provides the effective pore radius `r_eff`, the
porosity `epsilon`, the Darcy permeability `K` and the effective
liquid-saturated thermal conductivity `k_eff`.

References
----------
Faghri, *Heat Pipe Science and Technology*, Taylor & Francis, 1995.
Dunn & Reay, *Heat Pipes*, 5th ed., Elsevier, 2006.
"""
module Wick

export WickGeometry, screen_wick, sintered_wick, annular_wick,
       k_eff_maxwell, wick_mass

"""
    WickGeometry

Geometric/transport description of a capillary wick.
"""
struct WickGeometry
    kind::Symbol
    r_eff::Float64   # effective meniscus radius, m
    K::Float64       # Darcy permeability, m^2
    epsilon::Float64 # porosity, -
    t::Float64       # radial thickness, m
    k_s::Float64     # solid (wire/particle) conductivity, W/(m K)
end

"""
    screen_wick(d_wire, mesh_N; k_s = 15.0)

Square woven screen.  `mesh_N` is the number of wires per metre in each
direction (e.g. a 200 mesh screen has `N = 200/0.0254 = 7874 m^-1`).

* porosity: `epsilon = 1 - (pi/4) N d_wire`
* permeability (Blake-Kozeny): `K = d^2 epsilon^3 / (122 (1-epsilon)^2)`
* effective pore radius: `r_eff = 1/(2N)` (half the clear opening)
* thickness: `t = 2 d_wire`
"""
function screen_wick(d_wire::Real, mesh_N::Real; k_s::Real = 15.0)
    d = float(d_wire)
    N = float(mesh_N)
    eps = 1.0 - (pi / 4.0) * N * d
    eps = clamp(eps, 0.20, 0.95)
    K = d^2 * eps^3 / (122.0 * (1.0 - eps)^2)
    r_eff = 1.0 / (2.0 * N)
    return WickGeometry(:screen, r_eff, K, eps, 2.0 * d, float(k_s))
end

"""
    sintered_wick(d_particle; epsilon = 0.40, k_s = 15.0)

Sintered spherical powder wick.

* permeability (Kozeny-Carman): `K = d_p^2 epsilon^3 / (180 (1-epsilon)^2)`
* effective pore radius: `r_eff = 0.21 d_p`
* thickness defaults to `4 d_p`
"""
function sintered_wick(d_particle::Real; epsilon::Real = 0.40,
                       k_s::Real = 15.0, t::Union{Nothing,Real} = nothing)
    d = float(d_particle)
    eps = float(epsilon)
    K = d^2 * eps^3 / (180.0 * (1.0 - eps)^2)
    r_eff = 0.21 * d
    tt = t === nothing ? 4.0 * d : float(t)
    return WickGeometry(:sintered, r_eff, K, eps, tt, float(k_s))
end

"""
    annular_wick(t; r_pore = t/8, k_s = 15.0)

Annular gap of thickness `t`.  The liquid return path is the open annulus
(`K = t^2/12`), while the menisci are anchored on a fine screen or groove
surface whose pore radius `r_pore` controls the capillary pressure and the
entrainment limit.  This hybrid is representative of a self-venting or
screen-covered annular wick.
"""
function annular_wick(t::Real; r_pore::Real = float(t) / 8.0, k_s::Real = 15.0)
    tt = float(t)
    K = tt^2 / 12.0
    eps = 1.0
    return WickGeometry(:annular, float(r_pore), K, eps, tt, float(k_s))
end

"""
    k_eff_maxwell(k_l, k_s, epsilon)

Effective thermal conductivity of a liquid-saturated porous wick
(Maxwell / Faghri correlation):

    k_eff = k_l * ((k_l + k_s) - (1-eps)(k_l - k_s)) /
                  ((k_l + k_s) + (1-eps)(k_l - k_s))
"""
function k_eff_maxwell(k_l::Real, k_s::Real, epsilon::Real)
    num = (k_l + k_s) - (1.0 - epsilon) * (k_l - k_s)
    den = (k_l + k_s) + (1.0 - epsilon) * (k_l - k_s)
    return k_l * num / den
end

"""
    wick_mass(w, r_o, r_i, L, rho_screen)

Approximate dry wick mass in kg for a screen/sintered wick of thickness `t`
and solid density `rho_screen` (kg/m^3).
"""
function wick_mass(w::WickGeometry, r_o::Real, r_i::Real, L::Real, rho_screen::Real)
    area = pi * (r_i^2 - (r_i - w.t)^2)
    return (1.0 - w.epsilon) * area * L * rho_screen
end

end # module
