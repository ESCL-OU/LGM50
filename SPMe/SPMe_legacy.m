%% Single Particle Model w/ Electrolyte & Temperature simulated with CasADi
%   Legacy snapshot before SEI-growth reintroduction in SPMe.m
%   Preserved for regression comparison.

classdef SPMe_legacy
    properties
        % Parameters
        p
        % State Vectors
        y
        x
        u
        x_next
        % System Input
        cur
        % System Internal States
        c_ss_n
        c_ss_p
        c_s_n
        c_s_p
        c_s_n_dot
        c_s_p_dot
        c_n
        c_p
        c_e
        c_e_dot
        c_ex
        eta_p
        eta_pl
        jn
        jp
        jn_e
        jp_e
        % System Outputs
        V
        SOC
        % System Internal Functions
        f_n
        f_p
        f_e
        h_n
        h_p
        h
        % System Integrator
        f_int_n
        f_int_p
        f_int_e
        f_int
        % Simulator
        f_sim
        % Miscellaneous functions
        f_cex
        f_ss
        f_potentials
    end
    methods
        function obj = SPMe_legacy(p, solver)
            import casadi.*

            % Parameters
            obj.p = p;
            % System Internal States
            obj.c_s_n = MX.sym('c_s_n', obj.p.Nr-1);
            obj.c_s_p = MX.sym('c_s_p', obj.p.Nr-1);
            obj.c_n = MX.sym('c_n', obj.p.Nr+1);
            obj.c_p = MX.sym('c_p', obj.p.Nr+1);
            obj.c_e = MX.sym('c_e', obj.p.Nx-3);
            obj.c_e_dot = MX.sym('c_e_dot', length(obj.c_e));
            % System Input
            obj.cur = MX.sym('cur', 1);
            obj.u = obj.cur;
            % Ionic molar flux
            obj.jn = obj.cur / (obj.p.Faraday*obj.p.a_s_n*obj.p.Area_n*obj.p.L_n);
            obj.jp = -obj.cur / (obj.p.Faraday*obj.p.a_s_p*obj.p.Area_p*obj.p.L_p);

            % System Output variables
            obj.V = MX.sym('V');
            obj.SOC = MX.sym('SOC');
            obj.eta_pl = MX.sym('eta_pl');

            % State Vectors
            obj.x = vertcat(obj.c_s_n, obj.c_s_p, obj.c_e);

            % System Dynamics
            [A_n,A_p,B_n,B_p,C_n,C_p,D_n,D_p] = spm_plant_obs_mats(p);
            obj.c_s_n_dot = A_n*obj.c_s_n + B_n*obj.jn;
            obj.c_ss_n = C_n*obj.c_s_n + D_n*obj.jn;
            obj.c_s_p_dot = A_p*obj.c_s_p + B_p*obj.jp;
            obj.c_ss_p = C_p*obj.c_s_p + D_p*obj.jp;

            [obj.c_e_dot, obj.c_ex] = ode_electrolyte(obj.c_e, obj.jn, obj.jp, obj.p);

            % System Output
            [obj.V, obj.eta_pl, eta_n, eta_p, Unref, Upref, V_noVCE, V_electrolyteCond, V_electrolytePolar] = ...
                                  nonlinearSPMeOutputVoltage(obj.p, ...
                                                             obj.c_ss_n, obj.c_ss_p, ...
                                                             obj.c_ex, ...
                                                             obj.cur);
            obj.c_n = vertcat(obj.c_s_n(1), obj.c_s_n, obj.c_ss_n);
            obj.c_p = vertcat(obj.c_s_p(1), obj.c_s_p, obj.c_ss_p);

            % SOC Volume
            cs_ave_n = calc_c_ave(obj.c_n);
            obj.SOC = (cs_ave_n - obj.p.c_n_soc0) / (obj.p.c_n_soc100 - obj.p.c_n_soc0);

            obj.y = vertcat(obj.V, obj.SOC, obj.eta_pl, obj.c_ss_n);

            % System Internal Functions
            obj.f_n = Function('f_n', {obj.c_s_n, obj.cur}, {obj.c_s_n_dot}, {'x', 'p'}, {'x_dot'});
            obj.f_p = Function('f_p', {obj.c_s_p, obj.cur}, {obj.c_s_p_dot}, {'x', 'p'}, {'x_dot'});
            obj.f_e = Function('f_e', {obj.c_e, obj.cur}, {obj.c_e_dot}, {'x', 'p'}, {'x_dot'});
            obj.h_n = Function('h_n', {obj.c_s_n, obj.cur}, {obj.c_ss_n}, {'x', 'p'}, {'y'});
            obj.h_p = Function('h_n', {obj.c_s_p, obj.cur}, {obj.c_ss_p}, {'x', 'p'}, {'y'});

            % Integrator
            ode_n = struct('x', obj.c_s_n, 'p', obj.cur, 'ode', obj.c_s_n_dot);
            ode_p = struct('x', obj.c_s_p, 'p', obj.cur, 'ode', obj.c_s_p_dot);
            ode_e = struct('x', obj.c_e, 'p', obj.cur, 'ode', obj.c_e_dot);
            if ismember(solver, ['rk', 'rk4', 'idas', 'cvodes', 'collocation'])
                obj.f_int_n = integrator('f_int_n', solver, ode_n, 0, obj.p.delta_t);
                obj.f_int_p = integrator('f_int_p', solver, ode_p, 0, obj.p.delta_t);
                obj.f_int_e = integrator('f_int_e', solver, ode_e, 0, obj.p.delta_t);
            else
                obj.f_int_n = integrator('f_int_n', 'rk', ode_n, 0, obj.p.delta_t);
                obj.f_int_p = integrator('f_int_p', 'rk', ode_p, 0, obj.p.delta_t);
                obj.f_int_e = integrator('f_int_e', 'rk', ode_e, 0, obj.p.delta_t);
            end
            x_n_next = obj.f_int_n('x0', obj.c_s_n, 'p', obj.cur);
            x_p_next = obj.f_int_p('x0', obj.c_s_p, 'p', obj.cur);
            x_e_next = obj.f_int_e('x0', obj.c_e, 'p', obj.cur);
            obj.x_next = vertcat(x_n_next.xf, x_p_next.xf, x_e_next.xf);
            % System Function
            obj.f_int = Function('f_int', {obj.x, obj.u}, {obj.x_next}, {'x0', 'p'}, {'x_next'});
            obj.f_sim = Function('F', {obj.x, obj.u}, {obj.x_next, obj.y}, {'xk', 'u'}, {'xk_1', 'y'});
            obj.h = Function('h', {obj.x, obj.u}, {obj.y}, {'x', 'p'}, {'y'});
            obj.f_ss = Function('ss', {obj.x, obj.u}, {obj.c_ss_n, obj.c_ss_p}, {'x', 'u'}, {'c_ss_n', 'c_ss_p'});
            obj.f_cex = Function('c_ex', {obj.x}, {obj.c_ex}, {'x'}, {'c_ex'});
            obj.f_potentials = Function('f_potentials', {obj.x, obj.u}, ...
                {eta_n, eta_p, Unref, Upref, V_noVCE, V_electrolyteCond, V_electrolytePolar}, ...
                {'x', 'u'}, ...
                {'eta_n', 'eta_p', 'Unref', 'Upref', 'V_noVCE', 'V_electrolyteCond', 'V_electrolytePolar'});
        end

        function [x, y, c_ss_n, c_ss_p, c_ex] = step(obj, x0, I)
            F_int = obj.f_int('x0', x0, 'p', I);
            x = full(F_int.x_next);
            H = obj.h('x', x, 'p', I);
            y = full(H.y);
            c_ss = obj.f_ss('x', x, 'u', I);
            c_ss_n = full(c_ss.c_ss_n);
            c_ss_p = full(c_ss.c_ss_p);
            F_cex = obj.f_cex('x', x);
            c_ex = full(F_cex.c_ex);
        end

        function [eta_n, eta_p, Unref, Upref, V_noVCE, V_electrolyteCond, V_electrolytePolar] = potentials(obj, x, I)
            F_potential = obj.f_potentials('x', x, 'u', I);
            eta_n = full(F_potential.eta_n);
            eta_p = full(F_potential.eta_p);
            Unref = full(F_potential.Unref);
            Upref = full(F_potential.Upref);
            V_noVCE = full(F_potential.V_noVCE);
            V_electrolyteCond = full(F_potential.V_electrolyteCond);
            V_electrolytePolar = full(F_potential.V_electrolytePolar);
        end
    end
end
