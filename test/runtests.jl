using HeatPipeCore
using Test

@testset "HeatPipeCore" begin

    @testset "Fluid properties" begin
        for f in (Na, K, Li)
            Ts = range(T_melt(f) + 50.0, min(T_boil(f), T_crit(f) - 200.0),
                       length = 25)
            ps = [p_sat(f, T) for T in Ts]
            @test all(diff(ps) .> 0)
            @test all(p_sat(f, T) > 0 for T in Ts)
            @test all(h_fg(f, T) > 0 for T in Ts)
            @test all(rho_v(f, T) > 0 for T in Ts)
            @test all(rho_l(f, T) > 0 for T in Ts)
            @test all(sigma(f, T) > 0 for T in Ts)
            @test all(mu_l(f, T) > 0 for T in Ts)
            @test all(mu_v(f, T) > 0 for T in Ts)
        end
        # Sodium saturation pressure at the normal boiling point.
        @test p_sat(Na, T_boil(Na)) ≈ 1.01e5 rtol = 0.05
        # Latent heat decreases with temperature.
        @test h_fg(Na, 900.0) > h_fg(Na, 1200.0)
    end

    @testset "Wick correlations" begin
        w = screen_wick(1.0e-4, 2000.0)
        @test 0.2 < w.epsilon < 0.95
        @test w.K > 0
        @test w.r_eff > 0
        s = sintered_wick(5.0e-5)
        @test s.K > 0 && s.r_eff > 0
        a = annular_wick(5.0e-4)
        @test a.K ≈ (5.0e-4)^2 / 12.0
    end

    @testset "Operational limits" begin
        hp = krusty_heatpipe()
        for T in (773.15, 873.15, 973.15, 1073.15)
            lim = krusty_limit(hp, T)
            @test lim.minimum > 0
            @test isfinite(lim.minimum)
        end
        # Limits increase with temperature over the operating range.
        for fn in (q_viscous, q_sonic, q_entrainment)
            @test fn(hp, 1073.15) > fn(hp, 873.15)
        end
        # A space heat pipe has no gravity-assisted flooding branch.
        isr = isr_heatpipe(Na)
        @test q_flooding(isr, 1100.0) == Inf
        # Capillary pressure balance is positive for a horizontal pipe.
        @test q_capillary(isr, 1100.0) > 0
    end

    @testset "KRUSTY benchmark" begin
        v = krusty_validation()
        @test v.pass
        @test v.max_err <= 0.15
        for (T, err) in zip(v.T_C, v.rel_err)
            if T >= 550.0
                @test abs(err) <= 0.15
            end
        end
        np = nominal_point_validation()
        @test np.pass
        @test np.dT_ec <= 20.0
        @test abs(np.T_v_C - 800.0) <= 30.0
    end

    @testset "Walker entrainment benchmark" begin
        w = walker_validation()
        @test w.pass
        @test w.controlling == :entrainment
    end

    @testset "Steady-state network" begin
        hp = krusty_heatpipe(T_sink = 1063.0)
        r = steady_state(hp, 310.0; n = 60)
        @test r.T_v > hp.T_sink
        @test r.R_eq > 0
        @test r.dT_ec < 30.0
        # Mesh independence
        r30 = steady_state(hp, 310.0; n = 30)
        r120 = steady_state(hp, 310.0; n = 120)
        @test abs(r30.T_v - r120.T_v) < 2.0
    end

    @testset "Transient response" begin
        hp = krusty_heatpipe(T_sink = 1000.0)
        tr = transient_response(hp, t -> 310.0, (0.0, 400.0);
                                n = 20, dt = 0.1)
        @test length(tr.t) > 2
        @test all(isfinite, tr.T_v)
        @test tr.T_v[end] > 1000.0
        ss = steady_state(hp, 310.0; n = 20)
        @test abs(tr.T_v[end] - ss.T_v) < 30.0
    end

    @testset "ISR Mark 2 design space" begin
        c = isr_core()
        @test isr_pin_power(c) ≈ 8196.72 atol = 1.0
        f = isr_feasibility()
        @test maximum(r.Q_max for r in f) < 2000.0   # 0.6 cm cannot carry 8.2 kW
        @test all(r.Q_max < isr_pin_power(c) for r in f)

        a = isr_options()
        @test 0.01 < a.A.D_required < 0.04
        @test a.A.fuel_loss_pct > 5.0
        @test a.B.Q_total >= 0.95 * c.P_th
        @test a.C.Q_total < 0.25 * c.P_th

        m = isr_mass_balance()
        @test m.cooling_mass_saving_B_kg > 0
        @test m.pump_power_W > 0
        @test m.optionB_total_kg > m.mark1_total_kg
        @test m.optionA_total_kg < m.mark1_total_kg
    end

end
