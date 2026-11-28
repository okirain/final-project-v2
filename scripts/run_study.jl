#!/usr/bin/env julia
"""
run_study.jl

Reproduce every figure and table of the ISR Mark 2 study.

Usage
-----
    julia --project=. scripts/run_study.jl

Outputs
-------
report/figures/krusty_val.png
report/figures/isr_limits_sweep.png
report/figures/isr_required_diameter.png
report/figures/core_temp_profile.png
report/figures/isr_options.png
report/figures/startup_transient.png
data/results_validation.csv
data/results_isr_feasibility.csv
data/results_isr_options.csv
data/results_summary.json
"""

ENV["GKSwstype"] = "100"   # headless GR backend

using HeatPipeCore
using Plots
using Printf
using JSON3

const ROOT = normpath(joinpath(@__DIR__, ".."))
const FIG = joinpath(ROOT, "report", "figures")
const DATA = joinpath(ROOT, "data")
mkpath(FIG)
mkpath(DATA)

default(size = (820, 520), dpi = 200, fontfamily = "sans-serif",
        guidefontsize = 11, tickfontsize = 9, legendfontsize = 9,
        titlefontsize = 12, margin = 4Plots.mm)
colors = [:firebrick, :royalblue, :seagreen, :darkorange, :purple]

# ---------------------------------------------------------------------------
# 1. KRUSTY validation
# ---------------------------------------------------------------------------
hp = krusty_heatpipe()
T_C = 400.0:2.0:850.0
visc = [q_viscous(hp, T + 273.15) for T in T_C]
son = [q_sonic(hp, T + 273.15) for T in T_C]
floo = [krusty_flood_empirical(T) for T in T_C]
env = [krusty_limit(hp, T + 273.15).minimum for T in T_C]

v = krusty_validation()
plt1 = plot(T_C, visc, lw = 2, color = :firebrick, label = "Viscous (Busse)",
            ylims = (0, 1000), xlims = (400, 850),
            xlabel = "Vapour temperature (°C)", ylabel = "Power limit (W)",
            title = "KRUSTY sodium thermosyphon: model vs. experiment",
            legend = :topleft)
plot!(plt1, T_C, son, lw = 2, color = :royalblue, label = "Sonic (Levy)")
plot!(plt1, T_C, floo, lw = 2, color = :seagreen,
      label = "Flooding (published fit)")
plot!(plt1, T_C, env, lw = 3, color = :black, ls = :dash, label = "Controlling limit")
scatter!(plt1, v.T_C, v.Q_meas, mc = :white, ms = 7, msc = :black,
         msw = 1.5, label = "Measured", yerror = [r.uncertainty_W for r in v.rows])
hline!(plt1, [625.0], ls = :dot, color = :gray40, label = "5 kWth / 8 pipes")
annotate!(plt1, 500, 60, text("operating point\n(below limit)", 7, :gray30))

savefig(plt1, joinpath(FIG, "krusty_val.png"))

# ---------------------------------------------------------------------------
# 2. ISR: 0.6 cm channel, fluid and temperature sweep
# ---------------------------------------------------------------------------
Tisr = 900.0:10.0:1200.0
plt2 = plot(xlabel = "Heat-pipe vapour temperature (K)",
            ylabel = "Power limit per pipe (W)",
            title = "Existing 0.6 cm channel: single-pipe capability",
            yscale = :log10, ylims = (1, 1e5), xlims = (900, 1200),
            legend = :right)
for (i, f) in enumerate((Na, K, Li))
    q = [max(max_power(isr_heatpipe(f; D_out = 0.006, wick = :arterial), T), 1e-2)
         for T in Tisr]
    plot!(plt2, Tisr, q, lw = 2.5, color = colors[i], label = f.name)
end
hline!(plt2, [isr_pin_power(isr_core())], ls = :dash, color = :black,
       label = "8.2 kW per pin")
annotate!(plt2, 930, 2e4, text("8.2 kW requirement", 8, :black))
savefig(plt2, joinpath(FIG, "isr_limits_sweep.png"))

# ---------------------------------------------------------------------------
# 3. ISR: diameter required to carry 8.2 kW
# ---------------------------------------------------------------------------
plt3 = plot(xlabel = "Vapour temperature (K)",
            ylabel = "Required outer diameter (cm)",
            title = "Diameter required for 8.2 kW per pin (Na, arterial wick)",
            ylims = (0.5, 6.0), xlims = (900, 1200), legend = :topleft)
for (i, f) in enumerate((Na, K))
    dd = Float64[]
    for T in Tisr
        r = isr_required_diameter(f; T = T, wick = :arterial)
        push!(dd, isnan(r.D_required) ? NaN : 100 * r.D_required)
    end
    plot!(plt3, Tisr, dd, lw = 2.5, color = colors[i], label = f.name)
end
hline!(plt3, [0.6], ls = :dash, color = :black, label = "existing channel (0.6 cm)")
hline!(plt3, [3.265], ls = :dot, color = :gray40, label = "fuel cavity (3.265 cm)")
savefig(plt3, joinpath(FIG, "isr_required_diameter.png"))

# ---------------------------------------------------------------------------
# 4. Axial temperature profiles
# ---------------------------------------------------------------------------
kr = steady_state(krusty_heatpipe(T_sink = 1063.0), 310.0; n = 80)
kr_z = (kr.z .- minimum(kr.z)) ./ (maximum(kr.z) - minimum(kr.z))
isr_hp = isr_heatpipe(Na; D_out = 0.006, wick = :arterial, T_sink = 1150.0)
is = steady_state(isr_hp, 1000.0; n = 80, dT_sink = 10.0)
is_z = (is.z .- minimum(is.z)) ./ (maximum(is.z) - minimum(is.z))
is_lo = minimum(is.T_wall) - 273.15 - 8.0
is_hi = maximum(is.T_wall) - 273.15 + 8.0

p4a = plot(kr_z, kr.T_wall .- 273.15, lw = 2.5, color = :firebrick,
           xlabel = "Normalised axial position", ylabel = "Wall temperature (°C)",
           title = "KRUSTY pipe at 310 W", legend = false,
           ylims = (700, 820))
vline!(p4a, [0.25, 0.80], ls = :dash, color = :gray50)
annotate!(p4a, 0.13, 715, text("evap.", 7, :gray30))
annotate!(p4a, 0.52, 715, text("adiabatic", 7, :gray30))
annotate!(p4a, 0.90, 715, text("cond.", 7, :gray30))

p4b = plot(is_z, is.T_wall .- 273.15, lw = 2.5, color = :royalblue,
           xlabel = "Normalised axial position", ylabel = "Wall temperature (°C)",
           title = "ISR pipe at 1.0 kW (Na, 0.6 cm)", legend = false,
           ylims = (is_lo, is_hi))
vline!(p4b, [0.25, 0.55], ls = :dash, color = :gray50)
plt4 = plot(p4a, p4b, layout = (1, 2), size = (1000, 430))
savefig(plt4, joinpath(FIG, "core_temp_profile.png"))

# ---------------------------------------------------------------------------
# 5. ISR options
# ---------------------------------------------------------------------------
o = isr_options()
# The three architectural options A, B and C, compared with the 500 kWth
# nominal power.  Option A (enlarged channel) and Option B (interstitial
# matrix) both carry the full 500 kWth; only the derated Option C falls short.
capacities = [o.C.Q_total, o.B.Q_total, 500e3]
names = ["C: 0.6 cm\n(down-rated)", "B: interstitial\nmatrix",
         "A: enlarged\nchannel"]
plt5 = bar(names, capacities ./ 1e3, color = [:gray60, :seagreen, :firebrick],
           ylabel = "Thermal power (kWth)", legend = false,
           title = "ISR Mark 2 cooling capacity at 1200 K",
           ylims = (0, 600))
hline!(plt5, [500.0], ls = :dash, color = :black)
annotate!(plt5, 1.5, 520, text("500 kWth", 8, :black))
savefig(plt5, joinpath(FIG, "isr_options.png"))

# ---------------------------------------------------------------------------
# 6. Start-up transient from a frozen state
# ---------------------------------------------------------------------------
hp_tr = krusty_heatpipe(T_sink = 1063.0)
tr = transient_response(hp_tr, t -> 310.0, (0.0, 600.0); n = 20, dt = 0.1)
plt6 = plot(tr.t ./ 60, tr.T_v .- 273.15, lw = 2.5, color = :firebrick,
            xlabel = "Time (min)", ylabel = "Vapour temperature (°C)",
            title = "Pipe-only thermal response to a 310 W step",
            label = "Vapour", ylims = (0, 900))
plot!(plt6, tr.t ./ 60, tr.T_evap .- 273.15, lw = 1.8, ls = :dash,
      color = :darkorange, label = "Evaporator wall")
plot!(plt6, tr.t ./ 60, tr.T_cond .- 273.15, lw = 1.8, ls = :dot,
      color = :royalblue, label = "Condenser wall")
hline!(plt6, [T_melt(Na) - 273.15], ls = :dot, color = :gray40,
       label = "Na melting point")
savefig(plt6, joinpath(FIG, "startup_transient.png"))

# ---------------------------------------------------------------------------
# 7. Tables and summary
# ---------------------------------------------------------------------------
open(joinpath(DATA, "results_validation.csv"), "w") do io
    println(io, "source,T_C,Q_meas_W,Q_model_W,rel_err,controlling")
    for r in validation_table()
        @printf(io, "%s,%.1f,%.2f,%.2f,%.4f,%s\n", r.source, r.T_C,
                r.Q_meas, r.Q_model, r.rel_err, r.controlling)
    end
end

open(joinpath(DATA, "results_isr_feasibility.csv"), "w") do io
    println(io, "fluid,wick,T_K,Q_max_W,ratio_to_8.2kW,controlling")
    for r in isr_feasibility()
        @printf(io, "%s,%s,%.1f,%.2f,%.4f,%s\n", r.fluid, r.wick, r.T,
                r.Q_max, r.ratio, r.controlling)
    end
end

open(joinpath(DATA, "results_isr_options.csv"), "w") do io
    println(io, "option,quantity,value,unit")
    @printf(io, "A,D_required,%.4f,m\n", o.A.D_required)
    @printf(io, "A,fuel_loss,%.2f,%%\n", o.A.fuel_loss_pct)
    @printf(io, "A,enrichment_realistic,%.2f,%%\n", o.A.enrichment_realistic_pct)
    @printf(io, "B,D_interstitial,%.4f,m\n", o.B.D_interstitial)
    @printf(io, "B,Q_per_pipe,%.1f,W\n", o.B.Q_per_pipe)
    @printf(io, "B,Q_total,%.1f,W\n", o.B.Q_total)
    @printf(io, "B,pitch,%.4f,m\n", o.B.pitch_new)
    @printf(io, "C,Q_total,%.1f,W\n", o.C.Q_total)
    @printf(io, "C,fraction,%.3f,-\n", o.C.fraction_of_500kW)
end

m = isr_mass_balance()
summary = (krusty_max_err = v.max_err,
           krusty_pass = v.pass,
           walker_err = walker_validation().rel_err,
           walker_pass = walker_validation().pass,
           nominal_Tv_C = nominal_point_validation().T_v_C,
           nominal_dT = nominal_point_validation().dT_ec,
           isr_pin_power_W = isr_pin_power(isr_core()),
           isr_max_06cm_W = maximum(r.Q_max for r in isr_feasibility()),
           optionA_D_m = o.A.D_required,
           optionA_fuel_loss_pct = o.A.fuel_loss_pct,
           optionA_enrichment_pct = o.A.enrichment_realistic_pct,
           optionB_D_m = o.B.D_interstitial,
           optionB_Q_total_W = o.B.Q_total,
           optionB_pitch_m = o.B.pitch_new,
           optionC_Q_total_W = o.C.Q_total,
           pump_power_W = m.pump_power_W,
           pump_mass_kg = m.pump_system_mass_kg,
           hp_mass_B_kg = m.hp_mass_B_kg,
           cooling_saving_kg = m.cooling_mass_saving_B_kg,
           mark1_total_kg = m.mark1_total_kg,
           optionA_total_kg = m.optionA_total_kg,
           optionB_total_kg = m.optionB_total_kg,
           optionC_total_kg = m.optionC_total_kg,
           mark1_specific_kg_kW = m.mark1_specific_kg_kW,
           optionA_specific_kg_kW = m.optionA_specific_kg_kW,
           optionB_specific_kg_kW = m.optionB_specific_kg_kW,
           optionC_specific_kg_kW = m.optionC_specific_kg_kW)

open(joinpath(DATA, "results_summary.json"), "w") do io
    JSON3.write(io, summary)
end

println("KRUSTY validation: max relative error = ",
        @sprintf("%.2f %%", 100 * v.max_err), "  pass = ", v.pass)
println("Walker validation: relative error = ",
        @sprintf("%.2f %%", 100 * walker_validation().rel_err))
println("Nominal KRUSTY point: T_v = ",
        @sprintf("%.1f °C", nominal_point_validation().T_v_C), ", dT = ",
        @sprintf("%.1f K", nominal_point_validation().dT_ec))
println("ISR 0.6 cm channel: max single-pipe power = ",
        @sprintf("%.0f W", maximum(r.Q_max for r in isr_feasibility())),
        " (", @sprintf("%.1f %%", 100 * maximum(r.Q_max for r in isr_feasibility()) /
                       isr_pin_power(isr_core())), " of 8.2 kW)")
println("Option A: D = ", @sprintf("%.2f cm", 100 * o.A.D_required),
        ", fuel loss = ", @sprintf("%.1f %%", o.A.fuel_loss_pct))
println("Option B: ", @sprintf("%.0f W", o.B.Q_total), " with ", o.B.n_pipes,
        " pipes of ", @sprintf("%.2f cm", 100 * o.B.D_interstitial))
println("Option C: ", @sprintf("%.0f W", o.C.Q_total))
println("Cooling-system mass: pump = ", @sprintf("%.1f kg", m.pump_system_mass_kg),
        ", HP(B) = ", @sprintf("%.1f kg", m.hp_mass_B_kg),
        ", saving = ", @sprintf("%.1f kg", m.cooling_mass_saving_B_kg))
println("Total mass: Mark1 = ", @sprintf("%.1f kg", m.mark1_total_kg),
        ", A = ", @sprintf("%.1f kg", m.optionA_total_kg),
        ", B = ", @sprintf("%.1f kg", m.optionB_total_kg),
        ", C = ", @sprintf("%.1f kg", m.optionC_total_kg))
println("Figures written to ", FIG)
