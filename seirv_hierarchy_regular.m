function seirv_hierarchy_regular
% Closed (monokinetic) SEIRV system with hierarchy, integrated in the
% regular variables (rho_J^i, M_J^i = rho_J^i * m_J^i).
%
% With these variables no quotients rho_J'/rho_J appear, so there is no need
% for an eps_mass threshold or for deciding what m_J is when rho_J = 0.
% The only division is m_S^i = M_S^i / rho_S^i, and rho_S^i > 0 throughout
% the simulation.
%
% Equations (alpha_T = alpha_V = 1, variance of the susceptibles = 0):
%   theta = (beta_T/4) * sum_j f_j (rho_I^j - M_I^j)       [rho_I (1 - m_I)]
%   c_i   = theta * (rho_S^i - M_S^i)                      [infection, level i]
%   v_i   = (beta_V/2) * (rho_S^i + M_S^i)                 [vaccination, level i]
%   rho_S' = -c - v,  rho_E' = c - sigma_E rho_E,  rho_I' = sigma_E rho_E - gamma_I rho_I,
%   rho_R' = gamma_I rho_I,  rho_V' = v
%   M_S' = -(c+v) m_S + soc,   M_E' = c m_S - sigma_E M_E + soc,
%   M_I' = sigma_E M_E - gamma_I M_I + soc,  M_R' = gamma_I M_I + soc,
%   M_V' = v m_S + soc
% with the social term, for each compartment J of level i,
%   soc = (eta/tau) * [ p (m rho_J - M_J) + (1-p) (B_i rho_J - F_i M_J) ],
%   m = sum_j f_j sum_J M_J^j,   B_i = sum_{j>=i} f_j sum_J M_J^j,
%   F_i = sum_{j>=i} f_j.
%
% Reproduces Table 1 of the manuscript (values at t = 400) for the four
% scenarios. Common parameters: sigma_E = 0.2, gamma_I = 0.1, eta = 0.05,
% tau = 1, h = (0, 0.5, 1), rho_S^i(0) = 0.98, rho_I^i(0) = 0.02.

%% 1. SCENARIOS
% fields: name, bT (beta_T), bV (beta_V), p, f (level proportions), mS0 (initial
% opinion of the susceptibles at each level)
sc(1) = struct('name','1 Pure hierarchy',           'bT',0.40,'bV',0.01,'p',0.0,'f',[0.2 0.5 0.3],'mS0',[-0.5 0 0.5]);
sc(2) = struct('name','2 Effective vaccination', 'bT',0.20,'bV',0.15,'p',0.5,'f',[0.2 0.5 0.3],'mS0',[0.2 0.6 0.9]);
sc(3) = struct('name','3 Hesitant top level',    'bT',0.40,'bV',0.01,'p',0.0,'f',[0.4 0.4 0.2],'mS0',[0.8 0 -0.9]);
sc(4) = struct('name','4 High hierarchy breaking','bT',0.35,'bV',0.05,'p',0.9,'f',[0.2 0.5 0.3],'mS0',[-0.8 0 0.8]);

plot_case = 1;         % scenario to plot
T_final = 400;         % final time (days)
dt_out = 0.1;          % output step

%% 2. INTEGRATION AND TABLE
fprintf('\n%-24s %8s %9s %8s %9s %9s %9s %9s\n', 'Scenario', 'Reff(0)', ...
        'max rhoI', 't*', 'rhoR(400)', 'rhoV(400)', 'Reff(50)', 'm_inf');
res = cell(1, numel(sc));
for k = 1:numel(sc)
    P = set_parameters(sc(k));
    res{k} = simulate(P, T_final, dt_out);
    S = res{k};
    if S.kpeak == 1, tstar = '---'; else, tstar = sprintf('%.1f', S.t(S.kpeak)); end
    fprintf('%-24s %8.4f %9.4f %8s %9.4f %9.4f %9.4f %9.3f\n', P.name, ...
            S.Reff(1), S.peak, tstar, S.rhoR(end), S.rhoV(end), ...
            interp1(S.t, S.Reff, 50), S.minf);
end
fprintf('\nConsensus weights pi and eigenvalues of (eta/tau) L:\n');
for k = 1:numel(sc)
    S = res{k};
    fprintf('%-24s pi = [%s]   eigenvalues = [%s]\n', sc(k).name, ...
            num2str(S.pi', '%.3f '), num2str(S.eigvals', '%.4f '));
end

%% 3. PLOTS OF THE CHOSEN SCENARIO
S = res{plot_case}; P = set_parameters(sc(plot_case)); N = P.N;
figure(1);
plot(S.t, S.rhoS, 'b-', S.t, S.rhoE, 'g-', S.t, S.rhoI, 'r-', ...
     S.t, S.rhoR, 'c-', S.t, S.rhoV, 'm-', 'LineWidth', 2);
xlabel('Time (days)'); ylabel('Population fraction');
legend('S','E','I','R','V','Location','best'); grid on;
title(['Compartments (whole population): ' P.name]);

figure(2); hold on;
for i = 1:N
    plot(S.t, S.mi(:,i), 'LineWidth', 2, 'DisplayName', sprintf('h = %.1f', P.h(i)));
end
xlabel('Time (days)'); ylabel('m^i(t)'); legend('Location','best'); grid on;
title('Mean opinion by level');

figure(3); hold on;
for i = 1:N
    plot(S.t, S.rhoI_level(:,i), 'LineWidth', 2, 'DisplayName', sprintf('h = %.1f', P.h(i)));
end
xlabel('Time (days)'); ylabel('\rho_I^i(t)'); legend('Location','best'); grid on;
title('Infectious by level (fraction conditional on the level)');

figure(4);
plot(S.t, S.Reff, 'k-', S.t, ones(size(S.t)), 'k:', 'LineWidth', 2);
xlabel('Time (days)'); ylabel('R_{eff}(t)'); grid on;
title('Effective reproduction number');
drawnow;
end

%% ------------------------------------------------------------------------
function P = set_parameters(e)
% Fills in the parameters common to all scenarios.
P = e;
P.sE = 0.2;  P.gI = 0.1;          % day^-1
P.eta = 0.05; P.tau = 1;          % interaction strength and social time scale
P.h = [0 0.5 1];
P.N = numel(P.f);
assert(abs(sum(P.f) - 1) < 1e-12, 'The proportions must add up to 1.');
% Initial conditions (per level, conditional on the level)
P.r0 = zeros(P.N, 5);
P.r0(:,1) = 0.98;                 % S
P.r0(:,3) = 0.02;                 % I
% Initial opinion of each compartment. By default m_J^i(0) = m_S^i(0).
% It can be changed (for example, m_I different from m_S): P.m0 is N x 5.
P.m0 = repmat(P.mS0(:), 1, 5);
end

%% ------------------------------------------------------------------------
function S = simulate(P, T, dt)
N = P.N; f = P.f(:);
M0 = P.r0 .* P.m0;                % M_J^i(0) = rho_J^i(0) m_J^i(0)
y0 = [P.r0(:); M0(:)];
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[t, y] = ode45(@(t, y) rhs(t, y, P), 0:dt:T, y0, opts);   % not stiff with tau = 1
% If tau is small, replace ode45 by ode15s.

nt = numel(t);
R = reshape(y(:, 1:5*N),       nt, N, 5);   % R(:,i,J) = rho_J^i
M = reshape(y(:, 5*N+1:end),   nt, N, 5);   % M(:,i,J) = rho_J^i m_J^i
S.t = t;
S.rhoS = R(:,:,1) * f;  S.rhoE = R(:,:,2) * f;  S.rhoI = R(:,:,3) * f;
S.rhoR = R(:,:,4) * f;  S.rhoV = R(:,:,5) * f;
S.rhoI_level = R(:,:,3);
S.mi = sum(M, 3) ./ sum(R, 3);                  % mean opinion of level i
S.Reff = P.bT / (4 * P.gI) * ((R(:,:,1) .* (1 - S.mi).^2) * f);
[S.peak, S.kpeak] = max(S.rhoI);

% Consensus: matrix L (rows sum to zero), weights pi and limit opinion
L = zeros(N);
for i = 1:N
    for j = 1:N
        if j ~= i, L(i,j) = f(j) * (P.p + (1 - P.p) * (j > i)); end
    end
    L(i,i) = -sum(L(i,:));
end
pi_ = null(L');  pi_ = pi_ / sum(pi_);
S.pi = pi_;
S.eigvals = sort(real(eig(P.eta / P.tau * L)), 'descend');
S.minf = pi_' * P.mS0(:);                       % m_inf = pi . Mbar(0)
end


function dy = rhs(~, y, P)
N = P.N; f = P.f(:);
Fi = flipud(cumsum(flipud(f)));                 % F_i = sum_{j>=i} f_j
r = reshape(y(1:5*N),     N, 5);
M = reshape(y(5*N+1:end), N, 5);
rS = r(:,1); rI = r(:,3);
MS = M(:,1); ME = M(:,2); MI = M(:,3);

mS = MS ./ max(rS, 1e-14);                      % the only division in the system
theta = P.bT / 4 * sum(f .* (rI - MI));         % (bT/4) rho_I (1 - m_I)
c = theta * (rS - MS);                          % flow S -> E
v = P.bV / 2 * (rS + MS);                       % flow S -> V

dr = [ -c - v, ...
        c - P.sE * r(:,2), ...
        P.sE * r(:,2) - P.gI * rI, ...
        P.gI * rI, ...
        v ];

tot = sum(M, 2);                                % rho_J m_J summed over J, per level
m   = sum(f .* tot);                            % global mean opinion
B   = flipud(cumsum(flipud(f .* tot)));         % F_i b(t,h_i)
soc = (P.eta / P.tau) * ( P.p * (m * r - M) ...
      + (1 - P.p) * (bsxfun(@times, B, r) - bsxfun(@times, Fi, M)) );

dM = [ -(c + v) .* mS, ...
        c .* mS - P.sE * ME, ...
        P.sE * ME - P.gI * MI, ...
        P.gI * MI, ...
        v .* mS ] + soc;

dy = [dr(:); dM(:)];
end
