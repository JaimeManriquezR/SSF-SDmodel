function results = run_simulation(obj, model, options)
arguments
    obj (1,1) SDfilter
    model (1,1) SDmodel
    options.initial_conditions = []
    options.inflow_concentrations = []
    options.time_run (1,1) {mustBeNumeric} = 1;
    options.time_step (1,1) = "adaptive";
    options.nframes (1,1) {mustBeNumeric} = 200;
    options.phi_clogging (1,1) {mustBeNumeric} = .6;
    options.adaptive_velocity_factor (1,1) {mustBeNumeric} = .5;
    options.adaptive_time_tolerance (1,1) {mustBeNumeric} = .05;
end

initial_conditions = options.initial_conditions;
inflow_concentrations = options.inflow_concentrations;
time_run = options.time_run;
time_step = options.time_step;
nframes = options.nframes;
phi_clogging = options.phi_clogging;
dt_tolerance = options.adaptive_time_tolerance;
vel_safety = options.adaptive_velocity_factor;

if isnumeric(inflow_concentrations)
    C_in = inflow_concentrations;
    inflow_concentrations = @(t) [C_in(:)].';
end

%================= I. FILTER PARAMETERS ====================%
%-------  GEOMETRY -------------------------%
Hbar = -obj.domain(1);
B = obj.domain(2);

dz = obj.mesh.dz;
cell_centers = obj.mesh.cell_centers;
cell_boundaries = [-Hbar; cell_centers(1:end-1) + dz/2; B];

delta_index = cell_centers >= -obj.delta;
sand_index = cell_centers >= 0;
delta_index(sand_index) = false;

porosity = max(obj.epsilon,min(1,(obj.epsilon - 1)/obj.delta*cell_centers + obj.epsilon));
porosity_boundaries = max(obj.epsilon,min(1,(obj.epsilon - 1)/obj.delta*cell_boundaries + obj.epsilon));

NC = length(cell_centers);
n0 = find(abs(cell_centers) < dz/2);

%------- LIGHT ATTENUATION --------%
eta_water = model.light.attenuation.water*(cell_centers + Hbar);
eta_sand = 0*cell_centers;
eta_sand(delta_index) = (1 - obj.epsilon)*(cell_centers(delta_index) + obj.delta).*(1 + (cell_centers(delta_index) - obj.delta)/(2*obj.delta));
eta_sand(sand_index) = (1 - obj.epsilon)*(obj.delta/2 + cell_centers(sand_index));
eta_sand = model.light.attenuation.sand*eta_sand;

% ------ INFLOW VELOCITY -------%
q = obj.inflow_velocity./porosity_boundaries;

%========================================================%

%================= II. MODEL PARAMETERS ====================%
kP = length(model.particles.Name);
kL = length(model.liquids.Name);
rhoP = model.particles.density;
rhoL = model.liquids.density;

% assumption 1: all components (of the same type) have the same density
rhoP = rhoP(1);
rhoL = rhoL(1);

S_P = model.ecological.stoichiometric_matrix{model.particles.Name,:};
S_L = model.ecological.stoichiometric_matrix{model.liquids.Name,:};

alfa_matrix = diag([model.particles.dispersivity; ...
                    model.liquids.dispersivity]);

K_HetGrowth = model.ecological.kinetic_parameters{model.liquids.Name,"Het_Growth"};
K_PhoGrowth = model.ecological.kinetic_parameters{model.liquids.Name,"Pho_Growth"};
K_Hyd = model.ecological.kinetic_parameters{"DOM","Hydrolysis"};

att_diag = diag(model.particles.attachment_rate);
trP_diag = diag(model.particles.transfer_rate);
trL_diag = diag(model.liquids.transfer_rate);

mu = model.ecological.kinetic_parameters{"Reaction Rate",:}.*...
    (model.ecological.kinetic_parameters{"Temperature Correcting Factor",:}.^(293 - obj.temperature));
mu = mu(:); % get as column vector

fdark = model.ecological.dark_respiration;

% BIOFILM POROSITY
beta = model.beta;
% OSMOSIS RATE
tau = model.tau;
%========================================================%


%=================== III. PRE-ALLOCATION ======================%
%----------- CFL parameters -------------------%
alpha_P = model.particles.dispersivity(1);
alpha_L = model.liquids.dispersivity(1);

K_CFL = [K_HetGrowth, K_PhoGrowth];
att_max = max(att_diag(:));
trP_max = max(trP_diag(:));
trL_max = max(trL_diag(:));

% use nan for min-functions later, inf is used for CFL computation
K_HetGrowth(isinf(K_HetGrowth)) = nan;
K_PhoGrowth(isinf(K_PhoGrowth)) = nan;

%----------- SOLVER A -------------------------%
Id = speye(n0);
row_ = [2:n0 1:(n0-1) 2:n0 1:(n0-1)];
col_ = [(1:n0-1) 1:(n0-1) 2:n0 2:n0];

% CENTERED DIFFERENCES
S_q = spdiags([1 -1 0].*ones(n0,1),-1:1,n0,n0);
S_q(1) = -1;
S_q(end) = 0;

S_q = obj.inflow_velocity*(S_q./porosity(1:n0))/2/dz;

% LAPLACIAN MATRIX
D = spdiags([1 -2 1].*ones(n0,1),-1:1,n0,n0);
D = D/dz^2;

% ---------- GENERAL VARIABLES --------------%
Xi_list = cell(3,1); % used for Xi function (HET/POM monod)
phib_highest = 0;
vbmax_highest = 0;

frames_time = zeros(nframes,1);
concentration_frame_b = zeros(NC,nframes,kP + kP + kL);
concentration_frame_w = zeros(NC,nframes,1);
concentration_frame_f = zeros(NC,nframes,kP + kL);
reaction_frame_det = zeros(NC,nframes,kP);
reaction_frame_att_e = zeros(NC,nframes,kP);
reaction_frame_att_f = zeros(NC,nframes,kP);
reaction_frame_tr_P = zeros(NC,nframes,kP);
reaction_frame_tr_L = zeros(NC,nframes,kL);
%reaction_frame_rwe = zeros(NC,nframes,1);
reaction_frame_eco_M = zeros(NC,nframes,kP);
reaction_frame_eco_Pe = zeros(NC,nframes,kP);
reaction_frame_eco_Le = zeros(NC,nframes,kL);
reaction_frame_eco_Pf = zeros(NC,nframes,kP);
reaction_frame_eco_Lf = zeros(NC,nframes,kL);
light_frame = zeros(NC,nframes,1);

snap_counter = 1;
%========================================================%

%=================== IV. INITIAL CONDITIONS ======================%
inorganic_inflow = inflow_concentrations(0);
inorganic_inflow([1 2]) = 0;

default_conditions = struct( ...
    "biofilm", zeros(NC,2*kP + kL), ...
    "flowing", repmat(inorganic_inflow, [NC 1]), ...
    "enclosed_water", 0*cell_centers, ...
    "time_start", 0, ...
    "biofilm_velocity", [[];q(2:n0);0*q(n0+1:end-1);[]]);

% check for defined initial conditions and set default if not set
for field = string(fieldnames(default_conditions).')
    if ~isfield(initial_conditions,field)
        initial_conditions.(field) = default_conditions.(field);
    else
        if isempty(initial_conditions.(field))
            initial_conditions.(field) = default_conditions.(field);
        end
    end
end

% assign initial conditions
Cb = initial_conditions.biofilm;
Cf = initial_conditions.flowing;
phiWe = initial_conditions.enclosed_water;
time_start = initial_conditions.time_start;
vb = initial_conditions.biofilm_velocity;

% compute initial vf based on initial vb and phib(Cb,phiWe)
phib = sum(Cb(:,1:kP),2)/rhoP + ...
    sum(Cb(:,kP+1:kP+kP),2)/rhoP + ...
    sum(Cb(:,kP+kP+1:kP+kP+kL),2)/rhoL + phiWe;
phib_bdy = .5*(phib(2:end) + phib(1:end-1));

vf = [q(1);
    (q(2:end-1) - vb.*phib_bdy)./(1 - phib_bdy);
    q(end)];
%=========================================================================%

%=================== V. OUTPUT RESULTS ======================%
time_snap = linspace(time_start,time_start + time_run,nframes);

results = SDresults(obj,model);
results.time_start = time_start;
results.simulation_data.time_final_intended = time_start + time_run;
results.simulation_data.phi_clogging = phi_clogging;
results.simulation_data.biofilm_velocity = [obj.inflow_velocity, 0];
results.simulation_data.time = time_start;
results.simulation_data.CFL.value = [];
results.flag = "OK";

results.frames.reactions = array2table(zeros(5,3), ...
    'VariableNames',{'attachment', 'transfer', 'ecological'}, ...
    'rowNames',{'M', 'Pe', 'Le', 'Pf', 'Lf'}, ...
    'DimensionNames',{'Subphase','Reaction'});
%=========================================================================%

%=============== VI. TIME INTEGRATION ====================================%
if isnumeric(time_step)
    adaptivity = "none";
    dt = time_step;
    results.simulation_data.time_step = dt;
else
    adaptivity = time_step;
    dt = Inf;
    results.simulation_data.time_step = "time adapted";
end
%=========================================================================%

t = time_start;
while t < time_start + time_run
    C_in = inflow_concentrations(t);

    % ============= I. MAIN: computing state variables ================%

    %% GLOBAL CONCENTRATIONS
    % MATRIX
    CM = Cb(:,1:kP);
    % ENCLOSED REGION
    CPe = Cb(:,kP+1:kP+kP);
    CLe = Cb(:,kP+kP+1:kP+kP+kL);
    % FLOWING REGION
    CPf = Cf(:,1:kP);
    CLf = Cf(:,kP+1:kP+kL);

    %% VOLUME FRACTIONS
    % MATRIX
    phiM = sum(CM,2)/rhoP;
    % P. ENCLOSED + L. ENCLOSED + WATER
    phie = sum(CPe,2)/rhoP ...
        + sum(CLe,2)/rhoL ...
        + phiWe;
    % MATRIX + ENCLOSED
    phib = phiM + phie;
    phib_bdy = .5*(phib(2:end) + phib(1:end-1));
    % FLOWING SUSPENSION
    phif = 1 - phib;

    %% LOCAL CONCENTRATIONS
    % MATRIX
    Xb = CM./phib; Xb(phib == 0,:) = 0;
    Sb = CLe./phib; Sb(phib == 0,:) = 0;
    % ENCLOSED REGION
    Xe = CPe./phie; Xe(phie == 0,:) = 0;
    Se = CLe./phie; Se(phie == 0,:) = 0;
    % FLOWING REGION
    Xf = CPf./phif;
    Sf = CLf./phif;

    %% ATTACHMENT RATES (likelihood based on available volume)
    % these quantities are computed here in order to use them for CFL computing and save time
    a_e = phiM./phib; a_e(isnan(a_e)) = 0;
    a_f = (1 - porosity) + porosity.*phiM;

    vf_cell = .5*(vf(2:end) + vf(1:end-1));
    det_vf = model.detachment(vf_cell);


    %============= CHECK FOR FILTER CLOGGING =========================%
    if max(phib) > phib_highest
        phib_highest = max(phib);
        if phib_highest > phi_clogging
            CELL = find(phib > phi_clogging);
            warning('Clogged! %.0f%% biofilm!\nT = %f\nCELL = %i',100*phib_highest,t,CELL - n0)

            results.flag = "CLOGGED";
            results.simulation_data.error.description = "Biofilm volume fraction has surpassed clogging value.";
            results.simulation_data.error.problem_cells = CELL;
            results.simulation_data.error.time = t;

            break;
        end
    end
    %================================================================%

    %----- INEFFICIENT LOOP ---------------%
    loop_counter = 1;
    for X = {Xb, Xe, Xf}
        x = X{:}; % cell array to double array
        zeroCounter = (x(:,1).*x(:,3) == 0);
        Xi = zeros(size(cell_centers));
        Xi(~zeroCounter) = x(~zeroCounter,3)./(x(~zeroCounter,3) + K_Hyd*x(~zeroCounter,1));

        Xi_list{loop_counter} = Xi;
        loop_counter = loop_counter + 1;
    end
    [Xi_b, Xi_e, Xi_f] = deal(Xi_list{:});
    %--------------------------------------%

    %===================== II. TIME ADAPTIVITY ===========================================%
    switch adaptivity
        case "adaptive"

            ws_t0 = max(abs(S_P(1,3))*mu(3),abs(S_P(2,4))*mu(4));

            vbmax = (1 + vel_safety) * max(abs(vb));
            vfmax = (1 + vel_safety) * max(abs(vf));

            phib_max = max(phib);
            phie_f_max = max(phie./phif);

            ae_max = max(a_e);
            af_max = max(a_f);

            Xb_r = permute(Xb(:,[1 2]),[3 2 1]);
            Xe_r = permute(Xe(:,[1 2]),[3 2 1]);
            Xf_r = permute(Xf(:,[1 2]),[3 2 1]);

            %L_b = -sum(S_L_neg.*(mu'.*Xb_T(:,:,:))./(Sb_T + K_CFL),2);
            L_b = -sum(S_L(:,1:2).*(mu(1:2)'.*Xb_r)./(permute(Sb,[2 3 1]) + K_CFL),2);
            L_e = -sum(S_L(:,1:2).*(mu(1:2)'.*Xe_r)./(permute(Se,[2 3 1]) + K_CFL),2);
            L_f = -sum(S_L(:,1:2).*(mu(1:2)'.*Xf_r)./(permute(Sf,[2 3 1]) + K_CFL),2);

            w_v = [2*vbmax*[1, 1, 1], ...
                2*vfmax*[1, 1]];
            w_a = [0, 0, 0, ...
                2*vfmax*alpha_P*(1 + 1/(1 - phib_max)),...
                2*vfmax*alpha_L*(1 + 1/(1 - phib_max))];
            w_b = [max(det_vf), ...
                att_max*ae_max + trP_max/beta, ...
                max(trL_max/beta,1/tau), ...
                att_max*af_max + trP_max/beta*phie_f_max, ...
                trL_max/beta*phie_f_max];

            w_s = [max(ws_t0, abs(S_P(3,5))*mu(5)*max(Xi_b)), ...
                max(ws_t0, abs(S_P(3,5))*mu(5)*max(Xi_e)), ...
                max(abs(L_b),[],'all') + max(abs(L_e),[],'all'), ...
                max(ws_t0, abs(S_P(3,5))*mu(5)*max(Xi_f)), ...
                max(abs(L_f),[],'all')];

            dt_CFL = .99/max(w_v/dz + w_a/dz^2 + w_b + w_s);

            dt = min(dt_CFL,(1 + dt_tolerance)*dt);
            results.simulation_data.CFL.value(end+1) = dt_CFL;
    end
    %================================================================%

    %================= III. MAIN: compute reaction terms ===============================================%

    %% REACTION TERMS
    % LIGHT IRRADIATION
    eta_P = model.light.attenuation.particles*cumsum(sum(CM + CPe + CPf,2))*dz;
    I = obj.light_irradiation(t)*exp(-(eta_water + eta_sand + eta_P));

    % ATTACHMENT RX
    att_e = a_e.*CPe*att_diag;
    att_f = a_f.*CPf*att_diag;

    % DETACHMENT RX
    det_M = det_vf.*CM;

    % TRANSFER RX
    trns_P = (phie/beta).*(Xf - Xe)*trP_diag;
    trns_L = (phie/beta).*(Sf - Se)*trL_diag;

    monod_HetGrowth_b = min(Sb./(Sb + K_HetGrowth(:).'),[],2);
    monod_PhoGrowth_b = min(Sb./(Sb + K_PhoGrowth(:).'),[],2);

    monod_HetGrowth_e = min(Se./(Se + K_HetGrowth(:).'),[],2);
    monod_PhoGrowth_e = min(Se./(Se + K_PhoGrowth(:).'),[],2);

    monod_HetGrowth_f = min(Sf./(Sf + K_HetGrowth(:).'),[],2);
    monod_PhoGrowth_f = min(Sf./(Sf + K_PhoGrowth(:).'),[],2);

    rb = phib.*[mu(1)*monod_HetGrowth_b.*Xb(:,1), ...
        mu(2)*max(fdark,I.*exp(1 - I)).*monod_PhoGrowth_b.*Xb(:,2),...
        mu(3)*Xb(:,1), ...
        mu(4)*Xb(:,2), ...
        mu(5)*Xi_b.*Xb(:,1)];

    re = phie.*[mu(1)*monod_HetGrowth_e.*Xe(:,1), ...
        mu(2)*max(fdark,I.*exp(1 - I)).*monod_PhoGrowth_e.*Xe(:,2),...
        mu(3)*Xe(:,1), ...
        mu(4)*Xe(:,2), ...
        mu(5)*Xi_e.*Xe(:,1)];

    rf = phif.*[mu(1)*monod_HetGrowth_f.*Xf(:,1), ...
        mu(2)*max(fdark,I.*exp(1 - I)).*monod_PhoGrowth_f.*Xf(:,2),...
        mu(3)*Xf(:,1), ...
        mu(4)*Xf(:,2), ...
        mu(5)*Xi_f.*Xf(:,1)];

    % ECOLOGICAL RX
    sM = rb*S_P.';
    sPe = re*S_P.';
    sLe = (rb + re)*S_L.';
    sPf = rf*S_P.';
    sLf = rf*S_L.';

    % FULL RX TERMS
    RM = sM + att_e + att_f - det_M;
    RPe = sPe - att_e + trns_P;
    RLe = sLe + trns_L;
    RPf = sPf - att_f + det_M - trns_P;
    RLf = sLf - trns_L;

    % REACTION VECTORS
    Rb = [RM, RPe, RLe];
    Rf = [RPf, RLf];
    RWe = (beta*phib - phie)/tau;

    %============ SAVING FRAMES ================================%
    if t >= time_snap(snap_counter)
        frames_time(snap_counter) = t;

        concentration_frame_b(:,snap_counter,:) = Cb;
        concentration_frame_w(:,snap_counter,:) = rhoL*phiWe;
        concentration_frame_f(:,snap_counter,:) = Cf;

        reaction_frame_det(:,snap_counter,:) = det_M;
        reaction_frame_att_e(:,snap_counter,:) = att_e;
        reaction_frame_att_f(:,snap_counter,:) = att_f;
        reaction_frame_tr_P(:,snap_counter,:) = trns_P;
        reaction_frame_tr_L(:,snap_counter,:) = trns_L;
        %reaction_frame_rwe(:,snap_counter,:) = rhoL*RWe;
        reaction_frame_eco_M(:,snap_counter,:) = sM;
        reaction_frame_eco_Pe(:,snap_counter,:) = sPe;
        reaction_frame_eco_Le(:,snap_counter,:) = sLe;
        reaction_frame_eco_Pf(:,snap_counter,:) = sPf;
        reaction_frame_eco_Lf(:,snap_counter,:) = sLf;

        light_frame(:,snap_counter,:) = I;
        snap_counter = snap_counter + 1;
    end
    %==========================================================%

    %=================== IV. SOLVER A: compute biofilm velocity ======================%
    Rbio = sum(RM,2)/rhoP + sum(RPe,2)/rhoP + sum(RLe,2)/rhoL + RWe;

    % reassigning variables
    un = phib(1:n0);
    un_bdy = .5*(un(2:end) + un(1:end-1));

    % estimating Dirichlet data
    phib1 = phib(n0+1) + dt*Rbio(n0+1);

    % computing LHS
    l_n = model.cahn_hilliard.zeta_0*model.cahn_hilliard.mobility(un_bdy);
    l_n = l_n.*porosity_boundaries(2:n0);
    D_lambda = sparse(row_(:),col_(:),kron([1;-1],[l_n./porosity(2:n0); -l_n./porosity(1:n0-1)]))/dz^2;

    G = [Id - dt*S_q, -dt*D_lambda;
        -0/2*Id + model.cahn_hilliard.kappa*D,    Id];

    % computing RHS
    DPSI = model.cahn_hilliard.gradient_potential(un);

    RHS = [un + dt*Rbio(1:n0);
        -0/2*un + DPSI];

    RHS([n0+1;end]) = RHS([n0+1;end]) + (-model.cahn_hilliard.kappa/dz^2)*[0;phib1];

    % solving linear system
    sol = G\RHS;
    %un1 = sol(1:n0);
    %un1_bdy = .5*(un1(2:end) + un1(1:end-1));

    % computing vb
    vb = [[];
        q(2:n0) - model.cahn_hilliard.zeta_0*(1 - un_bdy).*(sol(n0+1+1:end) - sol(n0+1:end-1))/dz;
        0*q(n0+1:end-1);
        []];
    vbmax_highest = max(max(abs(vb)),vbmax_highest);

    %====================== V. SOLVER B: compute concentrations =================================%
    vf = [q(1);
        (q(2:end-1) - vb.*phib_bdy)./(1 - phib_bdy);
        q(end)];

    %%%% FLUX COMPUTING %%%%
    % BIOFILM
    Fb = [
        0*Cb(1,:);
        Cb(1:end-1,:).*max(0,vb) + Cb(2:end,:).*min(0,vb);
        0*Cb(1,:)
        ];
    Fin_b = porosity_boundaries(1:end-1).*Fb(1:end-1,:);
    Fout_b = porosity_boundaries(2:end).*Fb(2:end,:);

    % ENCLOSED WATER COMPONENT
    FWe = [
        0*phiWe(1,:);
        phiWe(1:end-1,:).*max(0,vb) + phiWe(2:end,:).*min(0,vb);
        0*phiWe(1,:)
        ];
    Fin_We = porosity_boundaries(1:end-1).*FWe(1:end-1,:);
    Fout_We = porosity_boundaries(2:end).*FWe(2:end,:);

    % FLOWING SUSPENSION
    varphi = Cf./phif;
    dnj = [
        0*Cf(1,:);
        dz\(abs(vf(2:end-2)).*(1 - phib_bdy(1:end-1))).*...
        ((varphi(2:end-1,:) - varphi(1:end-2,:))*alfa_matrix);
        0*Cf(1,:);
        0*Cf(1,:)
        ];

    Ff = [
        q(1)*C_in;
        Cf(1:end-1,:).*max(0,vf(2:end-1)) + Cf(2:end,:).*min(0,vf(2:end-1));
        q(end)*Cf(end,:)
        ] - dnj;

    Fin_f = porosity_boundaries(1:end-1).*Ff(1:end-1,:);
    Fout_f = porosity_boundaries(2:end).*Ff(2:end,:);

    results.simulation_data.biofilm_velocity(end+1,:) = [max(abs(vb)), max(sum(Fb,2))];

    %======================= VI. MAIN: update cell values =================================%
    Cb = Cb + (dt/dz)*(Fin_b - Fout_b)./porosity + dt*Rb;
    phiWe = phiWe + (dt/dz)*(Fin_We - Fout_We)./porosity + dt*RWe;
    Cf = Cf + (dt/dz)*(Fin_f - Fout_f)./porosity + dt*Rf;

    %========== CHECK IF CONCENTRATIONS ARE NEGATIVE ============%
    if any(Cb(:) < 0) || any(isnan(Cb(:)))
        CELL = mod(find(Cb(:) < 0),size(Cb,1));
        fprintf('Negative concentration in biofilm. \nTIME = %e\n',t)
        fprintf("CELL = %i\n",CELL(1) - n0)
        fprintf("HEIGHT = %i\n",cell_centers(CELL(1)));
        results.flag = "BIOFILM";

        results.simulation_data.error.description = "Concentration in biofilm volume has reached negative values.";
        results.simulation_data.error.problem_cells = CELL;
        results.simulation_data.error.time = t;

        break;
    elseif any(Cf(:) < 0) || any(isnan(Cf(:)))
        CELL = mod(find(Cf(:) < 0),size(Cf,1));
        fprintf('Negative concentration in flowing suspension. \nT = %e\n',t)
        fprintf("CELL = %i\n",CELL(1) - n0)
        fprintf("z = %i\n",cell_centers(CELL(1)));
        results.flag = "FLOWING";
        results.simulation_data.error.description = "Concentration in flowing suspension volume has reached negative values.";
        results.simulation_data.error.problem_cells = CELL;
        results.simulation_data.error.time = t;

        break;
    end
    %===========================================================%

    % NEXT STEP
    t = t + dt;
    results.simulation_data.time(end+1) = t;
end
frames_time(snap_counter) = t;

concentration_frame_b(:,snap_counter,:) = Cb;
concentration_frame_w(:,snap_counter,:) = rhoL*phiWe;
concentration_frame_f(:,snap_counter,:) = Cf;

reaction_frame_det(:,snap_counter,:) = det_M;
reaction_frame_att_e(:,snap_counter,:) = att_e;
reaction_frame_att_f(:,snap_counter,:) = att_f;
reaction_frame_tr_P(:,snap_counter,:) = trns_P;
reaction_frame_tr_L(:,snap_counter,:) = trns_L;
%reaction_frame_rwe(:,snap_counter,:) = rhoL*RWe;
reaction_frame_eco_M(:,snap_counter,:) = sM;
reaction_frame_eco_Pe(:,snap_counter,:) = sPe;
reaction_frame_eco_Le(:,snap_counter,:) = sLe;
reaction_frame_eco_Pf(:,snap_counter,:) = sPf;
reaction_frame_eco_Lf(:,snap_counter,:) = sLf;

light_frame(:,snap_counter,:) = I;

results.frames.time = frames_time;
results.simulation_data.finalvelocity_biofilm = vb;
results.simulation_data.finalvelocity_flowing = vf;

concentrations_cell = cell(kP + kL + 1,3);
for j = 1:kP
    concentrations_cell{j,1} = concentration_frame_b(:,:,j);
    concentrations_cell{j,2} = concentration_frame_b(:,:,kP + j);
    concentrations_cell{j,3} = concentration_frame_f(:,:,j);
end

for j = 1:kL
    concentrations_cell{kP + j,1} = 0;
    concentrations_cell{kP + j,2} = concentration_frame_b(:,:,kP + kP + j);
    concentrations_cell{kP + j,3} = concentration_frame_f(:,:,kP + j);
end
concentrations_cell{kP + kL + 1,1} = 0;
concentrations_cell{kP + kL + 1,2} = concentration_frame_w;
concentrations_cell{kP + kL + 1,3} = Inf;
results.frames.concentrations = cell2table(concentrations_cell, ...
    'VariableNames',{'matrix', 'enclosed', 'flowing'}, ...
    'rowNames',[results.model.particles.Name.', results.model.liquids.Name.', 'WATER'], ...
    'DimensionNames',{'Component','Volume'});

cell_attachment = {
    reaction_frame_det;
    reaction_frame_att_e;
    0;
    reaction_frame_att_f;
    0};
cell_transfer = {
    0;
    reaction_frame_tr_P;
    reaction_frame_tr_L;
    -reaction_frame_tr_P;
    -reaction_frame_tr_L};
cell_ecological = {
    reaction_frame_eco_M;
    reaction_frame_eco_Pe;
    reaction_frame_eco_Le;
    reaction_frame_eco_Pf;
    reaction_frame_eco_Lf};
results.frames.reactions = cell2table({cell_attachment{:}; cell_transfer{:}; cell_ecological{:}}.', ...
    'VariableNames',{'attachment', 'transfer', 'ecological'}, ...
    'rowNames',{'M', 'Pe', 'Le', 'Pf', 'Lf'}, ...
    'DimensionNames',{'Subphase','Reaction'});

results.frames.light = light_frame;

results.simulation_data.inflow_concentrations = inflow_concentrations;
results.simulation_data.phi_highest.value = phib_highest;
results.simulation_data.phi_highest.height = cell_centers(find(phib == phib_highest,1));
results.simulation_data.CFL.worst_case = model.CFL_bound(dz,max(q),obj.temperature,vbmax_highest,phi_clogging);
results.simulation_data.light_irradiation = obj.light_irradiation;

results.time_final = t;
end