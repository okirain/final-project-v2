"""
    ThermalNetwork

One-dimensional lumped-parameter thermal resistance network for a cylindrical
heat pipe.

The model resolves the envelope temperature along the axis (`n` finite-volume
nodes) and couples every node radially to a single, nearly isothermal vapour
node.  Axial conduction in the envelope, radial conduction through the wall
and the liquid-saturated wick, condensation film resistance and an external
sink conductance are represented explicitly.  The vapour node is treated as
quasi-steady (its thermal capacity is negligible compared with the metal).

Governing equations
-------------------
For every wall node `i`:

    C_i dT_i/dt = sum_j G_ax,ij (T_j - T_i)
                  + G_rad,i (T_v - T_i)
                  + q_i - G_sink,i (T_i - T_sink)

and for the vapour node:

    sum_i G_rad,i (T_i - T_v) = 0

At steady state the vapour balance is linear in `T_v` because the wall system
is linear, so the coupled problem is solved by a small fixed-point iteration
on the temperature-dependent properties.

The transient solver is a classical explicit RK4 integration of the wall
energy equations with the algebraic vapour temperature recomputed at every
stage.  This is sufficient for the lumped model and avoids a heavy ODE
dependency; the stable time step is reported by `stable_dt`.
"""
module ThermalNetwork

using ..FluidProperties
using ..Wick
using ..Limits

export ThermalNetworkModel, build_network, steady_state, transient_response,
       stable_dt, wall_capacity, axial_conductance

"""
    ThermalNetworkModel

Discretised thermal network.  All conductances are in W/K and capacities in
J/K.
"""
struct ThermalNetworkModel
    hp::HeatPipe
    n::Int
    z::Vector{Float64}        # node centre coordinates, m
    dz::Vector{Float64}       # node lengths, m
    region::Vector{Symbol}
    q::Vector{Float64}        # source per node, W
    G_sink::Vector{Float64}   # external sink conductance per node, W/K
    G_ax::Vector{Float64}     # axial conductance between i and i+1, W/K
    C::Vector{Float64}        # wall + wick heat capacity per node, J/K
    G_sink_total::Float64
end

"""
    build_network(hp; n = 60, Q = 0.0, dT_sink = 8.0, h_film = 5000.0)

Discretise `hp` into `n` equal axial finite volumes.  The source `Q` is
distributed uniformly over the evaporator and the external sink conductance
is sized so that the condenser wall sits `dT_sink` kelvin below the vapour at
the nominal power `Q`.
"""
function build_network(hp::HeatPipe; n::Int = 60, Q::Real = 0.0,
                       dT_sink::Real = 8.0, h_film::Real = 5000.0)
    L = total_length(hp)
    n = max(n, 6)
    dz = fill(L / n, n)
    z = [(i - 0.5) * L / n for i in 1:n]

    region = Vector{Symbol}(undef, n)
    for i in 1:n
        if z[i] <= hp.L_e
            region[i] = :evaporator
        elseif z[i] <= hp.L_e + hp.L_a
            region[i] = :adiabatic
        else
            region[i] = :condenser
        end
    end

    q = zeros(n)
    nevap = count(==(:evaporator), region)
    if nevap > 0 && Q > 0
        q[region .== :evaporator] .= float(Q) / nevap
    end

    ncond = count(==(:condenser), region)
    G_sink = zeros(n)
    G_sink_total = 0.0
    if ncond > 0 && Q > 0
        G_sink_total = float(Q) / max(float(dT_sink), 0.1)
        G_sink[region .== :condenser] .= G_sink_total / ncond
    end

    A_wall = pi * (hp.r_o^2 - hp.r_i^2)
    A_wick = pi * (hp.r_i^2 - hp.r_v^2)
    C = (hp.rho_wall * hp.cp_wall * A_wall .+
         0.5 * rho_l(hp.fluid, max(hp.T_sink, hp.fluid.T_melt)) *
         cp_l(hp.fluid, max(hp.T_sink, hp.fluid.T_melt)) * A_wick) .* dz

    # Harmonic mean of the two half-node conductances.
    G_ax = zeros(max(n - 1, 0))
    for i in 1:(n - 1)
        R1 = (dz[i] / 2.0) / (hp.k_wall * A_wall)
        R2 = (dz[i + 1] / 2.0) / (hp.k_wall * A_wall)
        G_ax[i] = 1.0 / (R1 + R2)
    end

    return ThermalNetworkModel(hp, n, z, dz, region, q, G_sink, G_ax, C,
                               G_sink_total)
end

"""Wall + wick heat capacity per node, J/K."""
wall_capacity(net::ThermalNetworkModel) = net.C

"""Axial envelope conductance vector, W/K."""
axial_conductance(net::ThermalNetworkModel) = net.G_ax

# ---------------------------------------------------------------------------
# Radial coupling
# ---------------------------------------------------------------------------

"""
    radial_conductance(net, i, T_v; h_film = 5000.0, frozen_factor = true)

Radial wall-to-vapour conductance of node `i`, W/K.
"""
function radial_conductance(net::ThermalNetworkModel, i::Int, T_v::Real;
                            h_film::Real = 5000.0, frozen_factor::Bool = true)
    hp = net.hp
    dz = net.dz[i]
    R_wall = log(hp.r_o / hp.r_i) / (2.0 * pi * hp.k_wall * dz)

    if net.region[i] == :condenser
        # Thin condensation film; h_film is a representative alkali-metal
        # film coefficient.  The film resistance is small relative to the wall.
        R_film = 1.0 / (float(h_film) * 2.0 * pi * hp.r_v * dz)
        G = 1.0 / (R_wall + R_film)
    else
        k_e = k_eff_maxwell(k_l(hp.fluid, T_v), hp.wick.k_s, hp.wick.epsilon)
        R_wick = log(hp.r_i / hp.r_v) / (2.0 * pi * k_e * dz)
        G = 1.0 / (R_wall + R_wick)
    end

    if frozen_factor
        # The working fluid cannot transport heat below its melting point;
        # switch the radial coupling on smoothly across a 20 K window.
        Tm = hp.fluid.T_melt
        G *= 0.02 + 0.98 / (1.0 + exp(-(T_v - Tm - 15.0) / 8.0))
    end
    return G
end

"""
    conductances(net, T_v) -> Vector{Float64}

All radial conductances evaluated at vapour temperature `T_v`.
"""
function conductances(net::ThermalNetworkModel, T_v::Real)
    return [radial_conductance(net, i, T_v) for i in 1:net.n]
end

# ---------------------------------------------------------------------------
# Tridiagonal solver
# ---------------------------------------------------------------------------

"""
    thomas!(a, b, c, d)

In-place Thomas algorithm for a tridiagonal system with sub-diagonal `a`
(length n-1), diagonal `b`, super-diagonal `c` (length n-1) and right-hand
side `d`.  Returns the solution vector.
"""
function thomas(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64},
                d::Vector{Float64})
    n = length(b)
    cp = zeros(n)
    dp = zeros(n)
    cp[1] = c[1] / b[1]
    dp[1] = d[1] / b[1]
    for i in 2:n
        m = b[i] - a[i - 1] * cp[i - 1]
        cp[i] = i < n ? c[i] / m : 0.0
        dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / m
    end
    x = zeros(n)
    x[n] = dp[n]
    for i in (n - 1):-1:1
        x[i] = dp[i] - cp[i] * x[i + 1]
    end
    return x
end

"""
    solve_wall(net, G_rad, T_v; ...) -> Vector{Float64}

Solve the steady wall energy balance for a given vapour temperature.
"""
function solve_wall(net::ThermalNetworkModel, G_rad::Vector{Float64},
                    T_v::Real, G_sink::Vector{Float64})
    n = net.n
    diag = zeros(n)
    lower = zeros(max(n - 1, 0))
    upper = zeros(max(n - 1, 0))
    rhs = net.q .+ G_rad .* T_v .+ G_sink .* net.hp.T_sink

    for i in 1:n
        diag[i] = G_rad[i] + G_sink[i]
        if i > 1
            diag[i] += net.G_ax[i - 1]
            lower[i - 1] = -net.G_ax[i - 1]
        end
        if i < n
            diag[i] += net.G_ax[i]
            upper[i] = -net.G_ax[i]
        end
    end
    return thomas(lower, diag, upper, rhs)
end

# ---------------------------------------------------------------------------
# Steady state
# ---------------------------------------------------------------------------

"""
    steady_state(hp, Q; n = 60, dT_sink = 8.0, tol = 1e-6, maxiter = 80)

Steady-state operating point of the heat pipe at throughput `Q`.

Returns a named tuple with the vapour temperature `T_v`, the axial wall
temperature profile `T_wall`, node coordinates `z`, the maximum
evaporator-to-condenser wall temperature difference `dT_ec`, the equivalent
thermal resistance `R_eq = (T_v - T_sink)/Q` and the limit set at `T_v`.
"""
function steady_state(hp::HeatPipe, Q::Real; n::Int = 60, dT_sink::Real = 8.0,
                      tol::Real = 1e-6, maxiter::Int = 200)
    net = build_network(hp; n = n, Q = Q, dT_sink = dT_sink)
    Q <= 0 && return (T_v = hp.T_sink, T_wall = fill(hp.T_sink, n), z = net.z,
                      dT_ec = 0.0, R_eq = 0.0, net = net,
                      limits = operating_limits(hp, hp.T_sink))

    # Vapour energy residual: net heat entering the vapour node.
    # It is monotone decreasing in T_v, so a bracketed bisection is robust
    # (the wall system is linear for fixed properties and is solved exactly).
    function residual(Tv)
        G = conductances(net, Tv)
        Tw = solve_wall(net, G, Tv, net.G_sink)
        return sum(G .* (Tw .- Tv))
    end

    lo = hp.fluid.T_melt + 1.0
    hi = min(0.98 * hp.fluid.T_crit, 2500.0)
    flo = residual(lo)
    fhi = residual(hi)
    it = 0
    while flo < 0.0 && it < 20
        lo = max(hp.fluid.T_melt * 0.5, lo - 50.0)
        flo = residual(lo)
        it += 1
    end
    it = 0
    while fhi > 0.0 && it < 20
        hi = min(hi + 100.0, 0.995 * hp.fluid.T_crit)
        fhi = residual(hi)
        it += 1
    end

    T_v = 0.5 * (lo + hi)
    for _ in 1:maxiter
        T_v = 0.5 * (lo + hi)
        fr = residual(T_v)
        if abs(fr) < 1e-6 || (hi - lo) < tol
            break
        end
        if sign(fr) == sign(flo)
            lo = T_v
            flo = fr
        else
            hi = T_v
        end
    end
    G_rad = conductances(net, T_v)
    T_wall = solve_wall(net, G_rad, T_v, net.G_sink)

    ev = net.region .== :evaporator
    cd = net.region .== :condenser
    dT_ec = (any(ev) && any(cd)) ? maximum(T_wall[ev]) - minimum(T_wall[cd]) : 0.0
    R_eq = Q > 0 ? (T_v - hp.T_sink) / Q : 0.0
    lims = operating_limits(hp, T_v)
    return (T_v = T_v, T_wall = T_wall, z = net.z, dT_ec = dT_ec,
            R_eq = R_eq, net = net, limits = lims)
end

# ---------------------------------------------------------------------------
# Transient
# ---------------------------------------------------------------------------

"""
    stable_dt(net) -> Float64

Conservative explicit-RK4 time step for the lumped wall network,
`dt = 0.25 min_i C_i / sum_j |G_ij|`.
"""
function stable_dt(net::ThermalNetworkModel)
    dt = Inf
    for i in 1:net.n
        Gsum = net.G_sink[i] + (i > 1 ? net.G_ax[i - 1] : 0.0) +
               (i < net.n ? net.G_ax[i] : 0.0)
        # Radial conductance at a representative high temperature.
        Gsum += radial_conductance(net, i, 900.0; frozen_factor = false)
        dt = min(dt, 0.25 * net.C[i] / max(Gsum, 1e-12))
    end
    return dt
end

"""
    transient_response(hp, Qfun, tspan; n = 40, dt = nothing, T0 = 300.0,
                       dT_sink = 8.0, store = 20)

Integrate the lumped network from a uniform initial temperature `T0`.

`Qfun(t)` returns the instantaneous source power in W (heat input to the
evaporator).  The external sink is held at `hp.T_sink`.  Returns a named
tuple `(t, T_v, T_evap, T_cond, T_wall)` where `T_evap`/`T_cond` are the
mean evaporator and condenser wall temperatures and `T_wall` is the final
profile.
"""
function transient_response(hp::HeatPipe, Qfun, tspan::Tuple{Real,Real};
                            n::Int = 40, dt::Union{Nothing,Real} = nothing,
                            T0::Real = 300.0, dT_sink::Real = 8.0,
                            store::Int = 20)
    Qmax = maximum(float(Qfun(tt)) for tt in range(tspan[1], tspan[2], length = 51))
    net = build_network(hp; n = n, Q = max(Qmax, 1.0), dT_sink = dT_sink)
    # Rebuild the source with the actual Qfun later; the network geometry is
    # what matters here.
    h = dt === nothing ? stable_dt(net) : float(dt)
    h = clamp(h, 0.02, 10.0)

    T = fill(float(T0), n)
    T_v = float(T0)
    G0 = conductances(net, T_v)
    T_v = sum(G0 .* T) / max(sum(G0), 1e-12)

    t = float(tspan[1])
    tend = float(tspan[2])
    nsteps = max(1, ceil(Int, (tend - t) / h))
    h = (tend - t) / nsteps
    nsave = max(1, ceil(Int, nsteps / store))

    ts = Float64[t]
    Tvs = Float64[T_v]
    Tev = Float64[mean_region(T, net, :evaporator)]
    Tcd = Float64[mean_region(T, net, :condenser)]

    function deriv(Tstate, tv, Q)
        G_rad = conductances(net, tv)
        dT = similar(Tstate)
        for i in 1:n
            axial = 0.0
            if i > 1
                axial += net.G_ax[i - 1] * (Tstate[i - 1] - Tstate[i])
            end
            if i < n
                axial += net.G_ax[i] * (Tstate[i + 1] - Tstate[i])
            end
            dT[i] = (axial + G_rad[i] * (tv - Tstate[i]) +
                     net.q[i] * (Q / max(Qmax, 1e-12)) -
                     net.G_sink[i] * (Tstate[i] - hp.T_sink)) / net.C[i]
        end
        return dT
    end

    function vapour_temp(Tstate)
        Tbar = sum(Tstate) / length(Tstate)
        G_rad = conductances(net, Tbar)
        s = sum(G_rad)
        return s > 1e-12 ? sum(G_rad .* Tstate) / s : Tstate[1]
    end

    for step in 1:nsteps
        Q = float(Qfun(t))
        k1 = deriv(T, T_v, Q)
        T2 = T .+ 0.5h .* k1
        tv2 = vapour_temp(T2)
        k2 = deriv(T2, tv2, float(Qfun(t + 0.5h)))
        T3 = T .+ 0.5h .* k2
        tv3 = vapour_temp(T3)
        k3 = deriv(T3, tv3, float(Qfun(t + 0.5h)))
        T4 = T .+ h .* k3
        tv4 = vapour_temp(T4)
        k4 = deriv(T4, tv4, float(Qfun(t + h)))
        T = T .+ (h / 6.0) .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        t += h
        T_v = vapour_temp(T)
        if step % nsave == 0 || step == nsteps
            push!(ts, t)
            push!(Tvs, T_v)
            push!(Tev, mean_region(T, net, :evaporator))
            push!(Tcd, mean_region(T, net, :condenser))
        end
    end
    return (t = ts, T_v = Tvs, T_evap = Tev, T_cond = Tcd, T_wall = T)
end

function mean_region(T::Vector{Float64}, net::ThermalNetworkModel, r::Symbol)
    idx = findall(==(r), net.region)
    return isempty(idx) ? NaN : sum(T[idx]) / length(idx)
end

end # module
