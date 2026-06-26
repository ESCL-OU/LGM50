%% Single Particle Model w/ Electrolyte & Temperature simulated with CasADi
%   Created Febuary 24, 2025 by Wonoo Choo
%   Called by sim_spmet.m

classdef SPMe
    properties
        % Parameters
        p
        degradation_model
        electrolyte_model
        use_lithium_plating
        use_anode_porosity_loss
        % State Vectors
        % x_n
        % x_p
        % u_n
        % u_p
        % x_dot_n
        % x_dot_p
        y
        % y_n
        % y_p
        % x_e
        % u_e
        % x_dot_e
        % x_n_next
        % x_p_next
        % x_e_next
        x
        u
        x_next
        x_dot
        % xk_1
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
        delta_sei
        delta_sei_dot
        delta_li
        delta_li_dot
        epsilon_e_n_dyn
        epsilon_e_n_dot
        c_ec_sei % added for Safari's SEI model
        c_ec_sei_dot % added for Safari's SEI model
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
        f_n % x_n_dot = f(x_n, u_n)
        f_p % x_p_dot = f(x_p, u_p)
        f_e % x_e_dot = f(x_e, u_e)
        h_n
        h_p
        h % y = h(x_n, u_n)
        % h_p % y = h(x_p, u_p)
        % System Integrator
        f_int_n
        f_int_p
        f_int_e
        f_int % f_int: (x, u) => (f_int)
        % Simulator
        f_sim % f_sim: (x, u) => (f_int, y)
        % Miscellaneous functions
        f_cex
        f_ss
        f_potentials
        f_optional_outputs
        optional_output_names
    end
    methods 
        function obj = SPMe(p, solver)
            import casadi.*
            
            % Parameters
            obj.p = p;
            % Degradation model selection at runtime:
            %   p.degradation_model in {'none','ramadass','yang','safari'}
            if isfield(obj.p, 'degradation_model') && ~isempty(obj.p.degradation_model)
                obj.degradation_model = lower(char(obj.p.degradation_model));
            else
                obj.degradation_model = 'none';
            end
            if ~ismember(obj.degradation_model, {'none', 'ramadass', 'yang', 'safari'})
                error('p.degradation_model must be one of: none, ramadass, yang, safari.');
            end
            if isfield(obj.p, 'enable_lithium_plating') && ~isempty(obj.p.enable_lithium_plating)
                obj.use_lithium_plating = ~strcmp(obj.degradation_model, 'none') && logical(obj.p.enable_lithium_plating);
            else
                obj.use_lithium_plating = false;
            end
            if isfield(obj.p, 'enable_anode_porosity_loss') && ~isempty(obj.p.enable_anode_porosity_loss)
                obj.use_anode_porosity_loss = ~strcmp(obj.degradation_model, 'none') && logical(obj.p.enable_anode_porosity_loss);
            else
                obj.use_anode_porosity_loss = false;
            end
            if obj.use_lithium_plating
                missing_plating_params = {};
                if ~isfield(obj.p, 'i0_lpl') || isempty(obj.p.i0_lpl)
                    missing_plating_params{end+1} = 'i0_lpl'; %#ok<AGROW>
                end
                if ~isfield(obj.p, 'alph_lpl') || isempty(obj.p.alph_lpl)
                    missing_plating_params{end+1} = 'alph_lpl'; %#ok<AGROW>
                end
                if ~isfield(obj.p, 'M_li') || isempty(obj.p.M_li)
                    missing_plating_params{end+1} = 'M_li'; %#ok<AGROW>
                end
                if ~isfield(obj.p, 'rho_li') || isempty(obj.p.rho_li)
                    missing_plating_params{end+1} = 'rho_li'; %#ok<AGROW>
                end
                if ~isempty(missing_plating_params)
                    error('Li plating enabled but missing parameter(s): %s.', strjoin(missing_plating_params, ', '));
                end
            end
            if isfield(obj.p, 'electrolyte_model') && ~isempty(obj.p.electrolyte_model)
                obj.electrolyte_model = lower(char(obj.p.electrolyte_model));
            else
                obj.electrolyte_model = 'fdm';
            end
            if ~ismember(obj.electrolyte_model, {'fdm', 'fvm'})
                error('p.electrolyte_model must be one of: fdm, fvm.');
            end
            use_sei_state = ~strcmp(obj.degradation_model, 'none');
            use_safari_state = strcmp(obj.degradation_model, 'safari'); % added for Safari's SEI model
            use_plating_state = obj.use_lithium_plating;
            use_porosity_state = obj.use_anode_porosity_loss;
            if use_safari_state % added for Safari's SEI model
                % Number of interior nodes used to discretize EC diffusion through the SEI film.
                % State order for safari mode is [c_s_n; c_s_p; c_e; delta_sei; c_ec_sei].
                if ~isfield(obj.p, 'Nsei') || isempty(obj.p.Nsei)
                    obj.p.Nsei = 10;
                end
                if ~isfield(obj.p, 'epsilon_sei') || isempty(obj.p.epsilon_sei)
                    % Safari uses K_EC*c_EC,sln at the SEI/electrolyte boundary and notes
                    % K_EC can be interpreted as SEI porosity. Default 1 keeps Yang-style c0_ec.
                    obj.p.epsilon_sei = 0.05; % 1.0;
                end
                if ~isfield(obj.p, 'delta_sei_min') || isempty(obj.p.delta_sei_min)
                    obj.p.delta_sei_min = 1e-12;
                end
            else
                obj.p.Nsei = 0;
            end

            % System Internal States
            obj.c_s_n = MX.sym('c_s_n', obj.p.Nr-1);
            obj.c_s_p = MX.sym('c_s_p', obj.p.Nr-1);
            % obj.c_ss_n = MX.sym('c_ss_n', 1);
            % obj.c_ss_p = MX.sym('c_ss_p', 1);
            obj.c_n = MX.sym('c_n', obj.p.Nr+1);
            obj.c_p = MX.sym('c_p', obj.p.Nr+1);
            obj.c_e = MX.sym('c_e', obj.p.Nx-3);
            obj.c_e_dot = MX.sym('c_e_dot', length(obj.c_e));
            if use_sei_state
                obj.delta_sei = MX.sym('delta_sei', 1);
            end
            if use_plating_state
                obj.delta_li = MX.sym('delta_li', 1);
            end
            if use_porosity_state
                obj.epsilon_e_n_dyn = MX.sym('epsilon_e_n_dyn', 1);
            end
            if use_safari_state % added for Safari's SEI model
                obj.c_ec_sei = MX.sym('c_ec_sei', obj.p.Nsei);
                obj.c_ec_sei_dot = MX.sym('c_ec_sei_dot', obj.p.Nsei);
            end
            % System Input
            obj.cur = MX.sym('cur', 1);
            obj.u = obj.cur;
            % Total ionic molar flux (before SEI side-reaction split)
            jn_tot = obj.cur / (obj.p.Faraday*obj.p.a_s_n*obj.p.Area_n*obj.p.L_n);
            jp_tot = -obj.cur / (obj.p.Faraday*obj.p.a_s_p*obj.p.Area_p*obj.p.L_p);

            % System Output variables
            obj.V = MX.sym('V');
            obj.SOC = MX.sym('SOC');
            obj.eta_pl = MX.sym('eta_pl');

            % State Vectors
            obj.x = vertcat(obj.c_s_n, obj.c_s_p, obj.c_e);
            if use_sei_state
                obj.x = vertcat(obj.x, obj.delta_sei);
            end
            if use_plating_state
                obj.x = vertcat(obj.x, obj.delta_li);
            end
            if use_porosity_state
                obj.x = vertcat(obj.x, obj.epsilon_e_n_dyn);
            end
            if use_safari_state % added for Safari's SEI model
                obj.x = vertcat(obj.x, obj.c_ec_sei);
            end

            % System Dynamics
            [A_n,A_p,B_n,B_p,C_n,C_p,D_n,D_p] = spm_plant_obs_mats(p);
            obj.c_ss_n = C_n*obj.c_s_n + D_n*jn_tot;
            obj.c_ss_p = C_p*obj.c_s_p + D_p*jp_tot;
            % Remark: I am cheating slightly here. jn_tot should be jn, but doing so
            % imposes an algebraic equation. This forms a DAE. I am going to
            % approximate jn by jn_tot, which should be ok, since jn and jn_tot have a
            % error on the order of 0.001%.
            
            [M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p);
            
            %% ELECTROLYTE PHASE DYNAMICS
            
            % Compute Boundary Conditions
            c_e_bcs = C_ce * obj.c_e;
            
            % Separate and aggregate electrolyte concentration
            c_en = obj.c_e(1:(obj.p.Nxn-1));
            c_es = obj.c_e((obj.p.Nxn-1)+1:(obj.p.Nxn-1)+(obj.p.Nxs-1));
            c_ep = obj.c_e((obj.p.Nxn-1)+obj.p.Nxs : end);
            c_ex = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];
            obj.c_ex = c_ex;

            % Electrolyte concentration endpoints
            ce0n = c_ex(1);
            cens = c_ex(p.Nxn+1);
            cesp = c_ex(p.Nxn+p.Nxs+1);
            ce0p = c_ex(end);

            %% VOLTAGE OUTPUT

            % Average electrolyte concentrations
            cen_bar = sum(c_ex(1:p.Nxn+1)) / (p.Nxn+1);
            ces_bar = sum(c_ex((p.Nxn+1):(p.Nxn+p.Nxs+1))) / (p.Nxs+1);
            cep_bar = sum(c_ex((p.Nxn+p.Nxs+1):(p.Nxn+p.Nxs+p.Nxp+1))) / (p.Nxp+1);

            ce_bar = casadi.MX.ones(p.Nx, 1);
            ce_bar(1:p.Nxn) = cen_bar;
            ce_bar((p.Nxn+1):(p.Nxn+p.Nxs-1)) = ces_bar;
            ce_bar(p.Nxn+p.Nxs+1:end) = cep_bar;

            [i_0n,i_0p] = exch_cur_dens(p, obj.c_ss_n, obj.c_ss_p, ce_bar, p.T_amb);
            RTaF = (p.R*p.T_amb)/(p.alph*p.Faraday);
            Unref = refPotentialAnode(p, obj.c_ss_n/p.c_s_n_max);
            Upref = refPotentialCathode(p, obj.c_ss_p/p.c_s_p_max);
            
            if use_sei_state
                delta_sei = obj.delta_sei;
            else
                delta_sei = 0;
            end
            if use_plating_state
                delta_li = obj.delta_li;
            else
                delta_li = 0;
            end
            if use_porosity_state
                epsilon_e_n = obj.epsilon_e_n_dyn;
            else
                epsilon_e_n = obj.p.epsilon_e_n;
            end
            delta_film = delta_sei + delta_li;
            if use_sei_state
                % Algebraically equivalent to Yang's w_sei*delta_film/kappa_sei.
                R_film = delta_sei / obj.p.kappa_sei;
            else
                R_film = 0;
            end

            R_tot_n = p.R_f_n + R_film;
            R_tot_p = p.R_f_p;

            % Voltage
            eta_n = RTaF * asinh(obj.cur / (2*p.a_s_n*p.Area_n*p.L_n*i_0n(1)));
            eta_p = RTaF * asinh(-obj.cur / (2*p.a_s_p*p.Area_p*p.L_p*i_0p(end)));
            V_noVCE =  eta_p - eta_n + Upref - Unref ...
                - (R_tot_n/(p.a_s_n*p.L_n*p.Area_n) + R_tot_p/(p.a_s_p*p.L_p*p.Area_p))*obj.cur;
            
            % Overpotentials due to electrolyte subsystem
            kap_n = electrolyteCond(p, cen_bar, p.T_amb);
            kap_s = electrolyteCond(p, ces_bar, p.T_amb);
            kap_p = electrolyteCond(p, cep_bar, p.T_amb);
            
            % Bruggeman relationships
            kap_n_eff = kap_n * epsilon_e_n.^(p.brug);
            kap_s_eff = kap_s * p.epsilon_e_s.^(p.brug);
            kap_p_eff = kap_p * p.epsilon_e_p.^(p.brug);

            % Activity-coefficient correction toggle for electrolyte polarization.
            use_electrolyte_activity = false;
            if isfield(p, 'use_electrolyte_activity') && ~isempty(p.use_electrolyte_activity)
                use_electrolyte_activity = logical(p.use_electrolyte_activity);
            end
            if use_electrolyte_activity
                dfca_n = electrolyteAct(cen_bar, p.T_amb, p);
                dfca_s = electrolyteAct(ces_bar, p.T_amb, p);
                dfca_p = electrolyteAct(cep_bar, p.T_amb, p);
            else
                dfca_n = 0;
                dfca_s = 0;
                dfca_p = 0;
            end

            % Overpotential due to electrolyte conductivity
            V_electrolyteCond = (p.L_n/(2*kap_n_eff) + 2*p.L_s/(2*kap_s_eff) + p.L_p/(2*kap_p_eff))*obj.cur; ...
                    
            % Overpotential due to electrolyte polarization
            V_electrolytePolar = (2*p.R*p.T_amb)/(p.Faraday) * ...
                                    ((1+dfca_n) * (log(cens) - log(ce0n)) ...
                                    +(1+dfca_s) * (log(cesp) - log(cens)) ...
                                    +(1+dfca_p) * (log(ce0p) - log(cesp)));

            % Add 'em up!
            obj.V = V_noVCE + V_electrolyteCond + V_electrolytePolar;

            obj.eta_pl = eta_n + Unref;


            %% Aging Dynamics

            j_sei = 0;
            j_side_total = 0;
            j_plating = 0;
            obj.delta_li_dot = 0;
            obj.epsilon_e_n_dot = 0;
            delta_film_dot = 0;
            eta_n_side = 0;
            phi_se = 0;
            eta_s = 0;
            if isfield(obj.p, 'c0_ec') && ~isempty(obj.p.c0_ec)
                c_s_ec = obj.p.c0_ec; % added for Safari's SEI model
            else
                c_s_ec = 0;
            end
            safari_pe = 0; % added for Safari's SEI model
            safari_k = 0; % added for Safari's SEI model
            switch obj.degradation_model
                case 'none'
                    obj.delta_sei_dot = 0;
                    obj.jn = jn_tot;
                    obj.jp = jp_tot;

                case 'ramadass'
                    %   SEI Layer Growth model
                    %   Eqns Adopted from Ramadass et al (2004) [Univ of South Carolina]
                    %   "Development of First Principles Capacity Fade Model for Li-Ion Cells"
                    %   DOI: 10.1149/1.1634273
                    %   NOTE1: This model has NOT been validated experimentally by eCAL
                    %   NOTE2: We assume this submodel only applies to anode
                    eta_n_side = RTaF * asinh(obj.cur / (2*obj.p.a_s_n*obj.p.Area_n*obj.p.L_n*i_0n(1)));
                    Unref_side = refPotentialAnode(obj.p, obj.c_ss_n/obj.p.c_s_n_max);
                    R_side = obj.p.R_f_n + obj.delta_sei/obj.p.kappa_P;
                    phi_se = eta_n_side + Unref_side + obj.p.Faraday*R_side*jn_tot;
                    eta_s = phi_se - obj.p.Us - obj.p.Faraday*R_side*jn_tot;
                    j_sei = -obj.p.i0s/obj.p.Faraday * exp((-obj.p.alph*obj.p.Faraday)/(obj.p.R*obj.p.T_amb)*eta_s);
                    obj.delta_sei_dot = -obj.p.M_P/(2 * obj.p.rho_P) * j_sei;
                    j_side_total = j_sei;

                case 'yang'
                    % Yang SEI growth model; optional add-ons are applied below.
                    eta_n_side = RTaF * asinh(obj.cur / (2*obj.p.a_s_n*obj.p.Area_n*obj.p.L_n*i_0n(1)));
                    Unref_side = refPotentialAnode(obj.p, obj.c_ss_n/obj.p.c_s_n_max);
                    phi_se = eta_n_side + Unref_side + obj.p.Faraday*R_film*jn_tot;
                    eta_s = phi_se - obj.p.U_sei - obj.p.Faraday*R_film*jn_tot;
                    K = obj.p.k0_sei * exp((-obj.p.alph*obj.p.Faraday)/(obj.p.R*obj.p.T_amb)*eta_s);
                    j_sei_numel = -K * obj.p.c0_ec;
                    j_sei_denomi = 1 + K * delta_film / obj.p.D_ec;
                    j_sei = j_sei_numel/j_sei_denomi;
                    obj.delta_sei_dot = -obj.p.M_sei/(2 * obj.p.rho_sei) * j_sei;

                    c_s_ec = j_sei * delta_film / obj.p.D_ec + obj.p.c0_ec;
                    obj.p.i0s = obj.p.k0_sei * c_s_ec * obj.p.Faraday;

                    j_side_total = j_sei;

                case 'safari'
                    % Safari et al. (JES 2009) SEI model with dynamic EC diffusion
                    % through the SEI film. Here j_sei is the EC-consumption molar
                    % surface flux [mol/(m^2 s)], so delta_dot includes 1/2 stoichiometry.
                    eta_n_side = RTaF * asinh(obj.cur / (2*obj.p.a_s_n*obj.p.Area_n*obj.p.L_n*i_0n(1)));
                    Unref_side = refPotentialAnode(obj.p, obj.c_ss_n/obj.p.c_s_n_max);
                    R_film = obj.delta_sei/obj.p.kappa_sei;

                    % Safari Eq. 5 driving potential: phi_1 - R_SEI*i_t.
                    % In this SPMe convention eta_n_side + Unref_side approximates phi_1,
                    % and F*R_film*jn_tot is the interfacial-area ohmic film drop.
                    phi_se = eta_n_side + Unref_side + obj.p.Faraday*R_film*jn_tot;
                    eta_s = phi_se - obj.p.U_sei - obj.p.Faraday*R_film*jn_tot;  % added U_sei term (which was omitted in Safari2009) to match with Yang                     
                    K = obj.p.k0_sei * exp((-obj.p.alph*obj.p.Faraday)/(obj.p.R*obj.p.T_amb)*eta_s);
                    safari_k = K;

                    % Moving-boundary growth factor: delta_dot = gamma*K*c_s_ec.
                    gamma_sei = obj.p.M_sei/(2*obj.p.rho_sei);
                    nsei = obj.p.Nsei;
                    delta_eff = obj.delta_sei + obj.p.delta_sei_min;
                    dz = delta_eff/(nsei + 1);
                    c_outer = obj.p.epsilon_sei * obj.p.c0_ec;
                    c1 = obj.c_ec_sei(1);

                    % Inner boundary couples diffusion/advection flux to reaction:
                    % -D*(c1-cs)/dz + v*cs = -K*cs, with v=gamma*K*cs.
                    % Stable positive root of: gamma*K*cs^2 + (D/dz+K)*cs - (D/dz)*c1 = 0.
                    A_bc = obj.p.D_ec/dz;
                    B_bc = A_bc + K;
                    c_s_ec = (2*A_bc*c1)/(B_bc + sqrt(B_bc^2 + 4*gamma_sei*K*A_bc*c1));

                    j_sei = -K * c_s_ec;
                    obj.delta_sei_dot = -gamma_sei * j_sei;
                    v_sei = obj.delta_sei_dot;
                    safari_pe = v_sei*delta_eff/obj.p.D_ec;

                    % EC transport inside SEI: dc/dt = D*d2c/dz2 - v*dc/dz.
                    c_full = vertcat(c_s_ec, obj.c_ec_sei, c_outer);
                    c_dot_vec = MX.zeros(nsei, 1);
                    for kk = 1:nsei
                        c_left = c_full(kk);
                        c_mid = c_full(kk+1);
                        c_right = c_full(kk+2);
                        d2c = (c_right - 2*c_mid + c_left)/(dz^2);
                        dcdr = (c_right - c_left)/(2*dz);
                        c_dot_vec(kk) = obj.p.D_ec*d2c - v_sei*dcdr;
                    end
                    obj.c_ec_sei_dot = c_dot_vec;

                    j_side_total = j_sei;
            end
            if obj.use_lithium_plating
                eta_lpl = phi_se - obj.p.Faraday * R_film * jn_tot;
                j_plating = -obj.p.i0_lpl/obj.p.Faraday * ...
                    exp((-obj.p.alph_lpl*obj.p.Faraday)/(obj.p.R*obj.p.T_amb) * eta_lpl);
                obj.delta_li_dot = -obj.p.M_li/obj.p.rho_li * j_plating;
                j_side_total = j_side_total + j_plating;
            end
            delta_film_dot = obj.delta_sei_dot + obj.delta_li_dot;
            if obj.use_anode_porosity_loss
                obj.epsilon_e_n_dot = -obj.p.a_s_n * delta_film_dot;
            end
            if use_sei_state
                obj.jn = (abs(jn_tot) - abs(j_side_total)) .* sign(jn_tot);
                obj.jp = jp_tot;
                obj.jn_e = jn_tot;
                obj.jp_e = jp_tot;
            else
                obj.jn = jn_tot;
                obj.jp = jp_tot;
                obj.jn_e = jn_tot;
                obj.jp_e = jp_tot;
            end
            %% Solid phase dynamics 

            obj.c_s_n_dot = A_n*obj.c_s_n + B_n*obj.jn;
            obj.c_s_p_dot = A_p*obj.c_s_p + B_p*obj.jp;

            %% Electrolyte dynamics 

            % Compute Electrolyte Diffusion Coefficient and Derivative
            [D_en,dD_en] = electrolyteDe(p, c_en, p.T_amb);
            [D_es,dD_es] = electrolyteDe(p, c_es, p.T_amb);
            [D_ep,dD_ep] = electrolyteDe(p, c_ep, p.T_amb);
            
            % ADD BRUGGEMAN RELATION % Apr.22 2016 by Saehong Park
            D_en_eff = D_en .* epsilon_e_n.^(p.brug-1);
            dD_en_eff = dD_en .* epsilon_e_n.^(p.brug-1);
            D_es_eff = D_es .* p.epsilon_e_s.^(p.brug-1);
            dD_es_eff = dD_es .* p.epsilon_e_s.^(p.brug-1);
            D_ep_eff = D_ep .* p.epsilon_e_p.^(p.brug-1);
            dD_ep_eff = dD_ep .* p.epsilon_e_p.^(p.brug-1);
            
            if strcmp(obj.electrolyte_model, 'fdm')
                % Compute derivative
                c_en_dot = dD_en_eff.*(M1n*c_en + M2n*c_e_bcs(1:2)).^2 ...
                    + D_en_eff.*(M3n*c_en + M4n*c_e_bcs(1:2)) + ...
                    ((1-p.t_plus)*p.a_s_n/epsilon_e_n) * obj.jn_e * ones(p.Nxn-1, 1) - ...
                    (obj.epsilon_e_n_dot/epsilon_e_n) * c_en;
                c_es_dot = dD_es_eff.*(M1s*c_es + M2s*c_e_bcs(2:3)).^2 ...
                    + D_es_eff.*(M3s*c_es + M4s*c_e_bcs(2:3));
                c_ep_dot = dD_ep_eff.*(M1p*c_ep + M2p*c_e_bcs(3:4)).^2 ...
                    + D_ep_eff.*(M3p*c_ep + M4p*c_e_bcs(3:4)) + diag(M5p)*obj.jp_e;
                obj.c_e_dot = [c_en_dot; c_es_dot; c_ep_dot];
            else
                p_electrolyte = obj.p;
                p_electrolyte.epsilon_e_n = epsilon_e_n;
                [obj.c_e_dot, ~] = ode_electrolyte_conservative(obj.c_e, obj.jn_e, obj.jp_e, p_electrolyte);
                idx_n = 1:(p.Nxn-1);
                obj.c_e_dot(idx_n) = obj.c_e_dot(idx_n) - ...
                    (obj.epsilon_e_n_dot/epsilon_e_n) * c_en;
                c_en_dot = obj.c_e_dot(1:(p.Nxn-1));
                c_es_dot = obj.c_e_dot((p.Nxn-1)+1:(p.Nxn-1)+(p.Nxs-1));
                c_ep_dot = obj.c_e_dot((p.Nxn-1)+p.Nxs : end);
            end
            c_en_dot_bar = sum(c_en_dot) / (p.Nxn-1);
            c_es_dot_bar = sum(c_es_dot) / (p.Nxs-1);
            c_ep_dot_bar = sum(c_ep_dot) / (p.Nxp-1);
            
            %% System Output 
            obj.c_n = vertcat(obj.c_s_n(1), obj.c_s_n, obj.c_ss_n);
            obj.c_p = vertcat(obj.c_s_p(1), obj.c_s_p, obj.c_ss_p);

            % SOC Volume
            cs_ave_n = calc_c_ave(obj.c_n);
            cs_ave_p = calc_c_ave(obj.c_p);
            c_ss_n_dot = C_n*obj.c_s_n_dot;
            c_ss_p_dot = C_p*obj.c_s_p_dot;
            c_n_dot_for_avg = vertcat(obj.c_s_n_dot(1), obj.c_s_n_dot, c_ss_n_dot);
            c_p_dot_for_avg = vertcat(obj.c_s_p_dot(1), obj.c_s_p_dot, c_ss_p_dot);
            cs_ave_n_dot = calc_c_ave(c_n_dot_for_avg);
            cs_ave_p_dot = calc_c_ave(c_p_dot_for_avg);
            obj.SOC = (cs_ave_n - obj.p.c_n_soc0) / (obj.p.c_n_soc100 - obj.p.c_n_soc0);
            % 
            % SOC Bulk
            % r_vec = (0:p.delta_r_n:1)';
            % obj.SOC = 3/p.c_s_n_max * trapz(r_vec,r_vec.^2.*obj.c_n);
            % 
            % SOC mean
            % cs_ave_n = mean(obj.c_s_n);
            % obj.SOC = (cs_ave_n - obj.p.c_n_soc0) / (obj.p.c_n_soc100 - obj.p.c_n_soc0);

            obj.y = vertcat(obj.V, obj.SOC, obj.eta_pl, obj.c_ss_n);

            % System Internal Functions
            % (In SEI-enabled mode, these subsystem derivatives depend on full state).
            obj.f_n = Function('f_n', {obj.x, obj.cur}, {obj.c_s_n_dot}, {'x', 'p'}, {'x_dot'});
            obj.f_p = Function('f_p', {obj.x, obj.cur}, {obj.c_s_p_dot}, {'x', 'p'}, {'x_dot'});
            obj.f_e = Function('f_e', {obj.x, obj.cur}, {obj.c_e_dot}, {'x', 'p'}, {'x_dot'});
            obj.h_n = Function('h_n', {obj.c_s_n, obj.cur}, {obj.c_ss_n}, {'x', 'p'}, {'y'});
            obj.h_p = Function('h_n', {obj.c_s_p, obj.cur}, {obj.c_ss_p}, {'x', 'p'}, {'y'});
            % obj.h = Function('h', {obj.x_n, obj.u}, {obj.y_n}, {'x', 'p'}, {'y'});
            % obj.h_p = Function('h_p', {obj.x_p, obj.u}, {obj.y_p}, {'x', 'p'}, {'y'});

            % Integrator
            obj.x_dot = vertcat(obj.c_s_n_dot, obj.c_s_p_dot, obj.c_e_dot);
            if use_sei_state
                obj.x_dot = vertcat(obj.x_dot, obj.delta_sei_dot);
            end
            if use_plating_state
                obj.x_dot = vertcat(obj.x_dot, obj.delta_li_dot);
            end
            if use_porosity_state
                obj.x_dot = vertcat(obj.x_dot, obj.epsilon_e_n_dot);
            end
            if use_safari_state % added for Safari's SEI model
                obj.x_dot = vertcat(obj.x_dot, obj.c_ec_sei_dot);
            end
            ode_full = struct('x', obj.x, 'p', obj.cur, 'ode', obj.x_dot);
            if ismember(solver, ['rk', 'rk4', 'idas', 'cvodes', 'collocation'])
                F_int = integrator('f_int_full', solver, ode_full, 0, obj.p.delta_t);
            else
                F_int = integrator('f_int_full', 'rk', ode_full, 0, obj.p.delta_t);
            end
            X_next = F_int('x0', obj.x, 'p', obj.cur);
            obj.x_next = X_next.xf;
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
            obj.optional_output_names = { ...
                'V', 'SOC', 'eta_pl', 'c_ss_n', 'c_ss_p', ...
                'delta_sei', 'delta_sei_dot', 'jn', 'jp', ...
                'j_sei', 'j_side', 'j_plating', 'i0s', ...
                'delta_li', 'delta_li_dot', 'delta_film', 'delta_film_dot', ...
                'epsilon_e_n', 'epsilon_e_n_dot', 'R_film', ...
                'eta_n_side', 'phi_se', 'eta_s', ...
                'c_s_ec', 'safari_pe', 'safari_k', ... % added for Safari's SEI model
                'V_electrolytePolar', 'ce0n', 'cens', 'cesp', 'ce0p', ...
                'dfca_n', 'dfca_s', 'dfca_p', ...
                'cen_bar', 'ces_bar', 'cep_bar', ...
                'cs_ave_n', 'cs_ave_p', ...
                'c_en_dot_bar', 'c_es_dot_bar', 'c_ep_dot_bar', ...
                'cs_ave_n_dot', 'cs_ave_p_dot'};
            obj.f_optional_outputs = Function('f_optional_outputs', {obj.x, obj.u}, ...
                {vertcat(obj.V, obj.SOC, obj.eta_pl, obj.c_ss_n, obj.c_ss_p, ...
                         delta_sei, obj.delta_sei_dot, obj.jn, obj.jp, ...
                         j_sei, j_side_total, j_plating, obj.p.i0s, ...
                         delta_li, obj.delta_li_dot, delta_film, delta_film_dot, ...
                         epsilon_e_n, obj.epsilon_e_n_dot, R_film, ...
                         eta_n_side, phi_se, eta_s, ...
                         c_s_ec, safari_pe, safari_k, ... % added for Safari's SEI model
                         V_electrolytePolar, ce0n, cens, cesp, ce0p, ...
                         dfca_n, dfca_s, dfca_p, ...
                         cen_bar, ces_bar, cep_bar, ...
                         cs_ave_n, cs_ave_p, ...
                         c_en_dot_bar, c_es_dot_bar, c_ep_dot_bar, ...
                         cs_ave_n_dot, cs_ave_p_dot)}, ...
                {'x', 'u'}, {'out'});
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

        function out = optional_outputs(obj, x, I, names)
            if nargin < 4 || isempty(names)
                names = obj.optional_output_names;
            elseif ischar(names) || isstring(names)
                names = {char(names)};
            end
            [tf, idx] = ismember(names, obj.optional_output_names);
            if ~all(tf)
                bad = strjoin(names(~tf), ', ');
                error('Unknown SPMe optional output names: %s', bad);
            end
            F_opt = obj.f_optional_outputs('x', x, 'u', I);
            v = full(F_opt.out);
            out = struct();
            for i = 1:numel(names)
                out.(names{i}) = v(idx(i));
            end
        end

    end 
end