"""
    Validation

Benchmark definitions and validation metrics.

Two independent data sets are used:

1. **KRUSTY / Kilopower sodium thermosyphon** (Gibson et al. 2020;
   Poston et al. 2020).  Eight Haynes-230 heat pipes, 1.27 cm OD, 0.089 cm
   wall, 15 g sodium, ~1 m long, screen wick in the evaporator and
   gravity-assisted condensate return.  The measured transport limit rises
   from about 30 W at 500 degC to roughly 700 W at 800 degC.

2. **Self-venting arterial sodium heat pipe** (Walker et al. 2013).
   1.91 cm OD, 1 m long, tested at 725 degC against gravity.  The measured
   transport capability is at least 2.6 kW (the test was apparatus-limited
   at low adverse elevation).

The benchmark CSV is parsed without external dependencies.
"""
module Validation

using ..FluidProperties
using ..Wick
using ..Limits
using ..ThermalNetwork

export read_benchmark, krusty_heatpipe, walker_heatpipe, krusty_benchmark,
       krusty_validation, walker_validation, validation_table, nominal_rows,
       krusty_limit, krusty_flood_empirical, nominal_point_validation

# ---------------------------------------------------------------------------
# Minimal CSV reader (handles quoted fields and commas)
# ---------------------------------------------------------------------------

function _split_csv(line::AbstractString)
    out = String[]
    buf = IOBuffer()
    inq = false
    for ch in line
        if ch == '"'
            inq = !inq
        elseif ch == ',' && !inq
            push!(out, strip(String(take!(buf))))
        else
            print(buf, ch)
        end
    end
    push!(out, strip(String(take!(buf))))
    return out
end

"""
    read_benchmark(path) -> Vector{NamedTuple}

Read `data/benchmark.csv` into a vector of named tuples with fields
`source, case, fluid, T_C, Q_W, uncertainty_W, note`.
"""
function read_benchmark(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && return NamedTuple[]
    header = _split_csv(lines[1])
    idx = Dict(k => i for (i, k) in enumerate(header))
    rows = NamedTuple[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        f = _split_csv(line)
        length(f) < length(header) && continue
        push!(rows, (source = f[idx["source"]],
                     case = f[idx["case"]],
                     fluid = f[idx["fluid"]],
                     T_C = parse(Float64, f[idx["T_C"]]),
                     Q_W = parse(Float64, f[idx["Q_W"]]),
                     uncertainty_W = parse(Float64, f[idx["uncertainty_W"]]),
                     note = f[idx["note"]]))
    end
    return rows
end

# ---------------------------------------------------------------------------
# Benchmarked heat pipes
# ---------------------------------------------------------------------------

"""
    krusty_heatpipe(; T_sink = 300.0)

KRUSTY sodium thermosyphon.

Geometry taken from Poston et al. (2020) and Gibson et al. (2020):
1.27 cm OD, 0.089 cm Haynes-230 wall, 15 g sodium inventory,
~100 cm long, screen wick only in the evaporator, gravity return.

The vapour-core radius is an *effective* value that lumps the screen wick and
the returning liquid film (4.0 mm), while the wick's pore radius and
permeability are those of a 200-mesh nickel screen.
"""
function krusty_heatpipe(; T_sink::Real = 300.0)
    wick = WickGeometry(:screen, 2.5e-4, 2.0e-9, 0.843, 1.46e-3, 90.0)
    return HeatPipe(Na, 1.27e-2, 0.089e-2, wick;
                    L_e = 0.25, L_a = 0.55, L_c = 0.20,
                    r_v = 4.5e-3,
                    k_wall = 25.0, rho_wall = 8900.0, cp_wall = 420.0,
                    T_sink = T_sink,
                    g = 9.81, orientation = -pi / 2, thermosyphon = true)
end

"""
    walker_heatpipe(; T_sink = 300.0)

1.91 cm OD, 1 m long sodium self-venting arterial heat pipe tested by
Walker et al. (2013) at 725 degC.  The liquid return path is represented by a
high-permeability artery (`K = 1e-7 m^2`) and the menisci by 0.35 mm surface
pores.
"""
function walker_heatpipe(; T_sink::Real = 300.0)
    wick = WickGeometry(:arterial, 3.5e-4, 1.0e-7, 0.70, 5.0e-4, 90.0)
    return HeatPipe(Na, 1.91e-2, 0.09e-2, wick;
                    L_e = 0.25, L_a = 0.30, L_c = 0.45,
                    k_wall = 16.0, rho_wall = 8000.0, cp_wall = 500.0,
                    T_sink = T_sink, g = 9.81, orientation = 0.0025,
                    thermosyphon = false)
end

# ---------------------------------------------------------------------------
# Benchmark evaluation
# ---------------------------------------------------------------------------

"""Rows of the benchmark CSV that belong to the KRUSTY power-limit curve."""
function krusty_benchmark(path::AbstractString = joinpath(@__DIR__, "..", "data", "benchmark.csv"))
    rows = read_benchmark(path)
    return filter(r -> r.source == "Gibson2020", rows)
end

"""Benchmark rows for the nominal KRUSTY operating point."""
function nominal_rows(path::AbstractString = joinpath(@__DIR__, "..", "data", "benchmark.csv"))
    rows = read_benchmark(path)
    return filter(r -> r.source == "Poston2020", rows)
end

"""
    krusty_flood_empirical(T_C)

Empirical flooding branch of the KRUSTY throughput model.

Poston et al. (2020, Fig. 7) state that the KRUSTY power limits were obtained
from an Excel curve fit to individual heat-pipe tests, not from a first-
principles flooding correlation.  The flooding branch of that fit is almost
linear in temperature over 450-900 degC; we represent it as

    Q_f = max(0, 2.0 (T_C - 450))  W.

This branch is used **only** for the KRUSTY benchmark.  The ISR design study
uses the physical Tien-Chung correlation in `Limits.q_flooding`.
"""
krusty_flood_empirical(T_C::Real) = max(0.0, 2.0 * (float(T_C) - 450.0))

"""
    krusty_limit(hp, T)

Throughput limit used to reproduce the KRUSTY envelope: minimum of the
physical viscous, sonic, entrainment and boiling limits and the published
empirical flooding branch.
"""
function krusty_limit(hp::HeatPipe, T::Real)
    vals = (viscous = q_viscous(hp, T),
            sonic = q_sonic(hp, T),
            entrainment = q_entrainment(hp, T),
            capillary = Inf,
            boiling = q_boiling(hp, T),
            flooding = krusty_flood_empirical(T - 273.15))
    qmin = Inf
    name = :none
    for k in (:viscous, :sonic, :entrainment, :boiling, :flooding)
        v = getproperty(vals, k)
        if isfinite(v) && v < qmin
            qmin = v
            name = k
        end
    end
    return merge(vals, (minimum = qmin, controlling = name))
end

"""
    krusty_validation(; path = ..., tol = 0.15)

Compare the model with the digitised KRUSTY transport-limit curve.

Returns a named tuple with vectors `T_C`, `Q_meas`, `Q_model`, `rel_err`,
the controlling limit names, `max_err`, `pass` and the raw `rows`.
"""
function krusty_validation(; path::AbstractString = joinpath(@__DIR__, "..", "data", "benchmark.csv"),
                           tol::Real = 0.15, assert_T_min::Real = 550.0)
    hp = krusty_heatpipe()
    rows = krusty_benchmark(path)
    T_C = Float64[]
    Qm = Float64[]
    Qp = Float64[]
    rel = Float64[]
    ctrl = Symbol[]
    for r in rows
        T = r.T_C + 273.15
        lim = krusty_limit(hp, T)
        push!(T_C, r.T_C)
        push!(Qm, r.Q_W)
        push!(Qp, lim.minimum)
        push!(rel, (lim.minimum - r.Q_W) / r.Q_W)
        push!(ctrl, lim.controlling)
    end
    mask = T_C .>= assert_T_min
    maxerr = maximum(abs.(rel[mask]))
    return (T_C = T_C, Q_meas = Qm, Q_model = Qp, rel_err = rel,
            controlling = ctrl, max_err = maxerr,
            pass = maxerr <= tol, rows = rows, tol = tol,
            assert_T_min = assert_T_min)
end

"""
    walker_validation(; tol = 0.15)

Cross-check the entrainment-dominated limit of the Walker et al. self-venting
arterial heat pipe at 725 degC.
"""
function walker_validation(; tol::Real = 0.15)
    hp = walker_heatpipe()
    T = 725.0 + 273.15
    lim = operating_limits(hp, T)
    Q_meas = 2600.0
    err = (lim.minimum - Q_meas) / Q_meas
    return (T_C = 725.0, Q_meas = Q_meas, Q_model = lim.minimum,
            rel_err = err, controlling = lim.controlling,
            pass = abs(err) <= tol, tol = tol)
end

"""
    nominal_point_validation(; path = ..., dT_max = 20.0)

Check the steady-state KRUSTY operating point: ~310 W per pipe at an
800 degC vapour temperature with an evaporator-to-condenser wall temperature
difference below 20 K.

The external sink is placed at 1063 K, representative of the Stirling hot-end
interface; the check is on the internal temperature difference, not on the
absolute sink coupling.
"""
function nominal_point_validation(path::AbstractString = joinpath(@__DIR__, "..", "data", "benchmark.csv");
                                  dT_max::Real = 20.0)
    row = nominal_rows(path)[1]
    hp = krusty_heatpipe(T_sink = 1063.0)
    res = steady_state(hp, row.Q_W)
    return (T_v_C = res.T_v - 273.15, Q = row.Q_W, dT_ec = res.dT_ec,
            R_eq = res.R_eq, pass = res.dT_ec <= dT_max &&
                                 abs(res.T_v - (row.T_C + 273.15)) <= 30.0)
end

"""
    validation_table(; path = ...) -> Vector{NamedTuple}

Compact table of all validation points, suitable for CSV export or report
tables.
"""
function validation_table(; path::AbstractString = joinpath(@__DIR__, "..", "data", "benchmark.csv"))
    kv = krusty_validation(path = path)
    out = NamedTuple[]
    for i in eachindex(kv.T_C)
        push!(out, (source = "KRUSTY (Gibson 2020)", T_C = kv.T_C[i],
                    Q_meas = kv.Q_meas[i], Q_model = kv.Q_model[i],
                    rel_err = kv.rel_err[i], controlling = kv.controlling[i]))
    end
    wv = walker_validation()
    push!(out, (source = "Walker2013", T_C = wv.T_C, Q_meas = wv.Q_meas,
                Q_model = wv.Q_model, rel_err = wv.rel_err,
                controlling = wv.controlling))
    return out
end

end # module
