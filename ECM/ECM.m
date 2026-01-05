classdef ECM
    properties
        % Parameters
        p % Battery Parameters
        ocv_p % OCV Curve Polynomial coefficients
        ecm_p % ECM model R/C values
        n_order 
        % State Vectors
        x
        x_dot
        x_next
        % System Input
        I
        % System Outputs
        V
        % System Internal Function
        f
        h
        % System Integator
        f_int
        % Simulator
        f_sim
    end
    methods
        function obj = ECM(p, ocv_p, ecm_p, n_order, solver)
            import casadi.*

            % Parameters
            obj.p = p;
            obj.ocv_p = ocv_p;
            obj.ecm_p = ecm_p;
            obj.n_order = n_order;

            % State Vector
            obj.x = MX.sym('x', n_order+1);
            % System Input
            obj.I = MX.sym('I', 1);
            % System Output
            obj.V = MX.sym('V');

            % System Dynamics
            obj.x_dot = ode_ecm(obj.x, obj.I, obj.p, obj.ecm_p, obj.n_order);
            obj.f = Function('f', {obj.x, obj.I}, {obj.x_dot}, {'x', 'I'}, {'x_dot'});
            % Integrator
            ode_f = struct('x', obj.x, 'p', obj.I, 'ode', obj.x_dot);
            if ismember(solver, ['rk', 'rk4', 'idas', 'cvodes', 'collocation'])
                f_int = integrator('f_int_', solver, ode_f, 0, obj.p.delta_t); 
            else
                f_int = integrator('f_int_', 'rk', ode_f, 0, obj.p.delta_t); 
            end
            x_next = f_int('x0', obj.x, 'p', obj.I);
            obj.x_next = x_next.xf;
            obj.f_int = Function('f_int', {obj.x, obj.I}, {obj.x_next}, {'x0', 'I'}, {'x_next'});
            obj.V = ecm_voltage(obj.x_next, obj.I, obj.ecm_p, obj.ocv_p, obj.n_order);
            obj.h = Function('h', {obj.x, obj.I}, {obj.V}, {'x', 'I'}, {'y'});
            obj.f_sim = Function('F_sim', {obj.x, obj.I}, {obj.x_next, obj.V}, {'x0', 'I'}, {'x_next', 'V'});
        end

        function [x, V] = step(obj, x0, I)
            F_sim = obj.f_sim('x0', x0, 'I', I);
            x = full(F_sim.x_next);
            V = full(F_sim.V);
        end
    end
end

