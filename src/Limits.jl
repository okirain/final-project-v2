"""
    Limits

Analytical steady-state operational limits of a cylindrical heat pipe:

1. viscous limit (Busse),
2. sonic limit (Levy),
3. entrainment / flooding limit (Cotter-Weber, Tien-Chung for thermosyphons),
4. capillary limit (Young-Laplace + Darcy / Hagen-Poiseuille),
5. boiling limit (Chi's nucleate-boiling criterion).

Each limit is returned in watts as a function of the vapour temperature `T`.
The heat pipe is represented by the `HeatPipe` geometry defined here.

The module deliberately avoids CFD: the vapour core is treated as a
one-dimensional laminar duct, the wick as a Darcy porous medium, and the
liquid/vapour interface as a Young-Laplace meniscus.
"""
module Limits

using ..FluidProperties
using ..Wick

export HeatPipe, vapour_area, wick_area, effective_length, total_length,
       hydraulic_diameter,
       q_viscous, q_sonic, q_entrainment, q_capillary, q_boiling, q_flooding,
       operating_limits, max_power, controlling_limit, limit_names,
       vapour_reynolds, is_sonic_choked

"""
    HeatPipe

Geometric and material description of a cylindrical heat pipe.

Fields are in SI units.  `orientation` is the tilt of the pipe axis above the
horizontal (radians), positive when the evaporator is above the condenser
(adverse, gravity opposes the liquid return) and negative when the condenser
is above the evaporator (assisting); `thermosyphon = true` enables the
gravity-assisted flooding limit.
"""
struct HeatPipe
    fluid::Fluid
    r_o::Float64          # envelope outer radius, m
    r_i::Float64          # envelope inner radius, m
    r_v::Float64          # vapour-core radius, m
    L_e::Float64          # evaporator length, m
    L_a::Float64          # adiabatic length, m
    L_c::Float64          # condenser length, m
    wick::WickGeometry
    k_wall::Float64       # envelope thermal conductivity, W/(m K)
    rho_wall::Float64     # envelope density, kg/m^3
    cp_wall::Float64      # envelope specific heat, J/(kg K)
    T_sink::Float64       # external sink temperature, K
    C_flood::Float64      # flooding coefficient; NaN -> Tien-Chung correlation
    g::Float64            # gravitational acceleration, m/s^2
    orientation::Float64  # radians, positive when evaporator is above condenser
    thermosyphon::Bool
end

function HeatPipe(f::Fluid, D_out::Real, wall_t::Real, wick::WickGeometry;
                  L_e::Real, L_a::Real, L_c::Real,
                  k_wall::Real = 25.0, rho_wall::Real = 8900.0,
                  cp_wall::Real = 420.0, T_sink::Real = 300.0,
                  C_flood::Real = NaN, g::Real = 0.0,
                  orientation::Real = 0.0, thermosyphon::Bool = false,
                  r_v::Union{Nothing,Real} = nothing)
    r_o = float(D_out) / 2.0
    r_i = r_o - float(wall_t)
    rv = r_v === nothing ? r_i - wick.t : float(r_v)
    rv <= 0.0 && throw(ArgumentError("vapour radius must be positive"))
    rv >= r_i && throw(ArgumentError("vapour radius must be smaller than r_i"))
    return HeatPipe(f, r_o, r_i, rv, float(L_e), float(L_a), float(L_c),
                    wick, float(k_wall), float(rho_wall), float(cp_wall),
                    float(T_sink), float(C_flood), float(g),
                    float(orientation), thermosyphon)
end

"""Vapour-core cross-sectional area, m^2."""
vapour_area(hp::HeatPipe) = pi * hp.r_v^2

"""Wick cross-sectional area available for liquid return, m^2."""
wick_area(hp::HeatPipe) = pi * (hp.r_i^2 - hp.r_v^2)

"""Effective heat-pipe length, m: `L_a + (L_e + L_c)/2`."""
effective_length(hp::HeatPipe) = hp.L_a + 0.5 * (hp.L_e + hp.L_c)

"""Total envelope length, m."""
total_length(hp::HeatPipe) = hp.L_e + hp.L_a + hp.L_c

"""Hydraulic diameter of the vapour core, m."""
hydraulic_diameter(hp::HeatPipe) = 2.0 * hp.r_v

# ---------------------------------------------------------------------------
# Individual limits
# ---------------------------------------------------------------------------

"""
    q_viscous(hp, T)

Busse viscous limit, W:

    Q_v = pi r_v^4 h_fg rho_v p_v / (16 mu_v l_eff)

The limit arises when the pressure gradient required to accelerate the
vapour exceeds the saturation-pressure gradient available along the pipe.
"""
function q_viscous(hp::HeatPipe, T::Real)
    f = hp.fluid
    p = p_sat(f, T)
    h = h_fg(f, T)
    rv = hp.r_v
    return pi * rv^4 * h * rho_v(f, T) * p /
           (16.0 * mu_v(f, T) * effective_length(hp))
end

"""
    q_sonic(hp, T)

Levy sonic limit, W:

    Q_s = A_v rho_v h_fg sqrt(gamma R_g T / (2 (gamma + 1)))

The vapour is choked (Mach = 1) at the evaporator exit.
"""
function q_sonic(hp::HeatPipe, T::Real)
    f = hp.fluid
    Rg = R_specific(f)
    a = sqrt(f.gamma * Rg * T / (2.0 * (f.gamma + 1.0)))
    return vapour_area(hp) * rho_v(f, T) * h_fg(f, T) * a
end

"""
    q_entrainment(hp, T)

Cotter entrainment (Weber-number) limit, W:

    Q_e = A_v h_fg sqrt(sigma rho_v / (2 r_h,w))

Liquid droplets are torn from the wick when the vapour shear stress exceeds
the surface-tension restoring force.
"""
function q_entrainment(hp::HeatPipe, T::Real)
    f = hp.fluid
    r_hw = hp.wick.r_eff
    return vapour_area(hp) * h_fg(f, T) *
           sqrt(sigma(f, T) * rho_v(f, T) / (2.0 * r_hw))
end

"""
    q_capillary(hp, T; include_gravity = true)

Capillary limit, W, obtained from the Young-Laplace / Darcy pressure balance

    dp_c,max = 2 sigma / r_eff >= dp_l + dp_v + dp_g

with

    dp_l = mu_l mdot l_eff / (rho_l K A_w)
    dp_v = 8 mu_v mdot l_eff / (pi rho_v r_v^4)
    dp_g = rho_l g L_total sin(orientation)

and `mdot = Q / h_fg`.  `orientation` is positive when the evaporator is above
the condenser (adverse), so `dp_g > 0` opposes the flow; the available driving
pressure is `dp_avail = dp_c,max - dp_g`.
"""
function q_capillary(hp::HeatPipe, T::Real; include_gravity::Bool = true)
    f = hp.fluid
    dp_c = 2.0 * sigma(f, T) / hp.wick.r_eff

    # Positive orientation = evaporator above condenser, so the hydrostatic
    # head opposes the liquid return (adverse); negative orientation assists.
    dp_g = 0.0
    if include_gravity
        dp_g = rho_l(f, T) * hp.g * total_length(hp) * sin(hp.orientation)
    end
    dp_avail = dp_c - dp_g
    dp_avail <= 0.0 && return 0.0

    h = h_fg(f, T)
    cl = mu_l(f, T) * effective_length(hp) /
         (rho_l(f, T) * hp.wick.K * wick_area(hp) * h)
    cv = 8.0 * mu_v(f, T) * effective_length(hp) /
         (pi * rho_v(f, T) * hp.r_v^4 * h)

    return dp_avail / (cl + cv)
end

"""
    q_boiling(hp, T; r_b = 2.54e-7)

Chi nucleate-boiling limit, W:

    Q_b = 2 pi L_e k_eff T_v / ln(r_i/r_v)
          * (2 sigma / (h_fg rho_v)) * (1/r_b - 1/r_eff)

`r_b` is the critical nucleation-bubble radius (default 0.254 micrometre).
"""
function q_boiling(hp::HeatPipe, T::Real; r_b::Real = 2.54e-7)
    f = hp.fluid
    k_e = k_eff_maxwell(k_l(f, T), hp.wick.k_s, hp.wick.epsilon)
    dT_drive = 2.0 * sigma(f, T) / (h_fg(f, T) * rho_v(f, T)) *
               (1.0 / r_b - 1.0 / hp.wick.r_eff)
    dT_drive <= 0.0 && return 0.0
    geom = 2.0 * pi * hp.L_e * k_e * T / log(hp.r_i / hp.r_v)
    return geom * dT_drive
end

"""
    q_flooding(hp, T)

Counter-current flooding limit of a gravity-assisted thermosyphon
(Tien-Chung / Faghri), W:

    Q_f = C_w^2 A_v h_fg rho_v^(1/2) [g sigma (rho_l - rho_v)]^(1/4)

with `C_w^2 = 3.2 tanh^2(0.5 Bo^(1/4))` and
`Bo = 2 r_v sqrt(g (rho_l - rho_v) / sigma)`.

Returns `Inf` when gravity is absent or the pipe is not a thermosyphon.
"""
function q_flooding(hp::HeatPipe, T::Real)
    (hp.g <= 0.0 || !hp.thermosyphon) && return Inf
    f = hp.fluid
    rl = rho_l(f, T)
    rv = rho_v(f, T)
    s = sigma(f, T)
    drho = max(rl - rv, 0.1 * rl)
    Bo = 2.0 * hp.r_v * sqrt(hp.g * drho / s)
    C2 = isnan(hp.C_flood) ? 3.2 * tanh(0.5 * Bo^0.25)^2 : hp.C_flood
    # Flooding requires the condenser above the evaporator (assist), i.e. a
    # negative orientation under the adverse-positive convention.
    C2 *= max(-sin(hp.orientation), 0.0)
    return C2 * vapour_area(hp) * h_fg(f, T) * sqrt(rv) *
           (hp.g * s * drho)^0.25
end

"""
    operating_limits(hp, T) -> NamedTuple

All limits evaluated at `T`, in watts, plus the controlling limit.
"""
function operating_limits(hp::HeatPipe, T::Real)
    qv = q_viscous(hp, T)
    qs = q_sonic(hp, T)
    qe = q_entrainment(hp, T)
    # In a gravity-assisted thermosyphon the condensate returns as a film, so
    # the capillary pressure balance does not set the throughput.
    qc = hp.thermosyphon ? Inf : q_capillary(hp, T)
    qb = q_boiling(hp, T)
    qf = q_flooding(hp, T)
    vals = (viscous = qv, sonic = qs, entrainment = qe,
            capillary = qc, boiling = qb, flooding = qf)
    qmin, name = findmin_named(vals)
    return merge(vals, (minimum = qmin, controlling = name))
end

const _LIMIT_ORDER = (:viscous, :sonic, :entrainment, :capillary, :boiling, :flooding)

function findmin_named(nt::NamedTuple)
    best = Inf
    name = :none
    for k in _LIMIT_ORDER
        v = getproperty(nt, k)
        if isfinite(v) && v < best
            best = v
            name = k
        end
    end
    return (best, name)
end

"""Minimum (controlling) power limit in W at temperature `T`."""
max_power(hp::HeatPipe, T::Real) = operating_limits(hp, T).minimum

"""Symbol naming the controlling limit at temperature `T`."""
controlling_limit(hp::HeatPipe, T::Real) = operating_limits(hp, T).controlling

"""Human-readable names of the limit fields."""
limit_names() = _LIMIT_ORDER

"""
    vapour_reynolds(hp, T, Q)

Vapour-core Reynolds number at throughput `Q`, used to flag transition from
the laminar (Hagen-Poiseuille) to the turbulent (Blasius) regime.
"""
function vapour_reynolds(hp::HeatPipe, T::Real, Q::Real)
    f = hp.fluid
    mdot = Q / h_fg(f, T)
    return 2.0 * mdot / (pi * hp.r_v * mu_v(f, T))
end

"""`true` when the vapour exit velocity reaches Mach 1 at power `Q`."""
function is_sonic_choked(hp::HeatPipe, T::Real, Q::Real)
    return Q >= q_sonic(hp, T)
end

end # module
