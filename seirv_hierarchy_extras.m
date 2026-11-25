function seirv_hierarchy_extras
% Figures of the sensitivity section of the manuscript:
%   (1) heat maps in the (p, m_S^3(0)) plane              -> heat_p_m3.jpg
%   (2) heat maps in the (p, tau) plane                    -> heat_p_tau.jpg
%   (3) closed system versus Monte Carlo (Example 3)       -> mc_vs_closure.jpg
%
% It uses the same closed system, in regular variables (rho_J^i, M_J^i),
% as seirv_hierarchy_regular.m. The Monte Carlo simulates the Boltzmann
% model (B) with agents.
%
% Each block can be switched off with the flags below. Rough running times:
% heat maps, a few minutes; Monte Carlo, several minutes (3 modes x 4 seeds).

do_heat1 = true;
do_heat2 = true;
do_mc    = true;

%% Parameters of the controlled design (R0 = 1.5 with neutral opinions)
B = struct('bT',0.6,'bV',0.02,'p',0,'f',[0.4 0.4 0.2],'mS0',[0.8 0 0], ...
           'sE',0.2,'gI',0.1,'eta',0.05,'tau',1);
T_heat = 600;

%% (1) Map in (p, m3)
if do_heat1
    ps  = 0:0.1:1;            % 11 values
    m3s = -1:0.2:1;           % 11 values
    Rf = zeros(numel(m3s), numel(ps));   % recovered fraction at T_heat
    Pk = zeros(numel(m3s), numel(ps));   % peak of rho_I
    for a = 1:numel(m3s)
        for b = 1:numel(ps)
            P = B; P.p = ps(b); P.mS0 = [0.8 0 m3s(a)];
            S = integrate(P, T_heat, 1, false);
            Rf(a,b) = S.rhoR(end);  Pk(a,b) = max(S.rhoI);
        end
    end
    plot_map(ps, m3s, Rf, Pk, 'p', 'm_S^3(0)', 'heat_p_m3.jpg');
    fprintf('\nMap (p, m3):\n');
    fprintf('  rho_R(%d) with p = 0: from %.3f to %.3f;  with p = 1: from %.3f to %.3f\n', ...
            T_heat, min(Rf(:,1)), max(Rf(:,1)), min(Rf(:,end)), max(Rf(:,end)));
    fprintf('  cells with an interior peak (peak > 0.0201):\n');
    for a = 1:numel(m3s)
        k = find(Pk(a,:) > 0.0201);
        if ~isempty(k)
            fprintf('    m3 = %5.1f : p = %s\n', m3s(a), mat2str(ps(k), 2));
        end
    end
end

%% (2) Map in (p, tau), with m3 = -0.9
if do_heat2
    ps   = 0:0.1:1;
    taus = logspace(-2, 1, 10);
    Rf2 = zeros(numel(taus), numel(ps));
    Pk2 = zeros(numel(taus), numel(ps));
    for a = 1:numel(taus)
        for b = 1:numel(ps)
            P = B; P.p = ps(b); P.mS0 = [0.8 0 -0.9]; P.tau = taus(a);
            S = integrate(P, T_heat, 1, true);    % small tau: stiff system
            Rf2(a,b) = S.rhoR(end);  Pk2(a,b) = max(S.rhoI);
        end
    end
    plot_map(ps, log10(taus), Rf2, Pk2, 'p', 'log_{10} \tau', 'heat_p_tau.jpg');
    fprintf('\nMap (p, tau), m3 = -0.9:\n');
    fprintf('  p = 0: rho_R from %.3f (tau = %.2f) to %.3f (tau = %.0f);  peak from %.3f to %.3f\n', ...
            Rf2(1,1), taus(1), Rf2(end,1), taus(end), Pk2(1,1), Pk2(end,1));
    fprintf('  p = 1: rho_R between %.3f and %.3f\n', min(Rf2(:,end)), max(Rf2(:,end)));
end

%% (3) Closed system versus Monte Carlo (Example 3)
if do_mc
    P3 = struct('bT',0.4,'bV',0.01,'p',0,'f',[0.4 0.4 0.2],'mS0',[0.8 0 -0.9], ...
                'sE',0.2,'gI',0.1,'eta',0.05,'tau',1);
    T = 400; dt = 0.1; Na = 1e5; nseed = 4;
    Scl = integrate(P3, T, dt, false);
    t = (0:dt:T)';
    modes = {'dirac', 'unif', 'bimodal'};
    names = {'Dirac at each level', 'uniform, half-width (0.2, 1, 0.1)', 'bimodal middle level at +-1'};
    fprintf('\nMonte Carlo of the Boltzmann model (%d agents, dt = %.2f, %d seeds):\n', Na, dt, nseed);
    fprintf('  closed system: peak %.4f, t* = %.1f, rho_R(400) = %.4f\n', ...
            max(Scl.rhoI), t(find(Scl.rhoI == max(Scl.rhoI), 1)), Scl.rhoR(end));
    Imean_dirac = []; Isd_dirac = [];
    for m = 1:numel(modes)
        I = zeros(numel(t), nseed); Rr = zeros(1, nseed);
        for s = 1:nseed
            [I(:,s), Rr(s)] = montecarlo(P3, s, modes{m}, T, dt, Na);
        end
        pk = max(I, [], 1);
        ts = zeros(1, nseed);
        for s = 1:nseed, ts(s) = t(find(I(:,s) == pk(s), 1)); end
        fprintf('  %-34s peak %.4f +- %.4f,  t* = %.1f +- %.1f,  rho_R(400) = %.4f +- %.4f\n', ...
                names{m}, mean(pk), std(pk), mean(ts), std(ts), mean(Rr), std(Rr));
        if m == 1, Imean_dirac = mean(I, 2); Isd_dirac = std(I, 0, 2); end
    end
    figure(10); clf; hold on;
    fill([t; flipud(t)], [Imean_dirac - Isd_dirac; flipud(Imean_dirac + Isd_dirac)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, Scl.rhoI, 'k-', 'LineWidth', 2, 'DisplayName', 'closed system');
    plot(t, Imean_dirac, 'r--', 'LineWidth', 2, 'DisplayName', 'Monte Carlo (mean of 4 runs)');
    xlabel('Time (days)'); ylabel('\rho_I(t)'); legend('Location','best'); grid on;
    print(gcf, '-djpeg', '-r130', 'mc_vs_closure.jpg');
end
end

%% ------------------------------------------------------------------------
function plot_map(xs, ys, Z1, Z2, xl, yl, filename)
figure;
subplot(1,2,1); imagesc([xs(1) xs(end)], [ys(1) ys(end)], Z1);
set(gca,'YDir','normal'); axis tight; colorbar;
xlabel(xl); ylabel(yl); title('recovered fraction \rho_R(600)');
subplot(1,2,2); imagesc([xs(1) xs(end)], [ys(1) ys(end)], Z2);
set(gca,'YDir','normal'); axis tight; colorbar;
xlabel(xl); ylabel(yl); title('epidemic peak max_t \rho_I');
print(gcf, '-djpeg', '-r130', filename);
end

%% ------------------------------------------------------------------------
function S = integrate(P, T, dt, stiff)
% Closed system in regular variables. Returns the global rho_I and rho_R.
N = numel(P.f); f = P.f(:);
r0 = zeros(N, 5); r0(:,1) = 0.98; r0(:,3) = 0.02;
M0 = r0 .* repmat(P.mS0(:), 1, 5);
y0 = [r0(:); M0(:)];
if stiff
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [t, y] = ode15s(@(t, y) rhs(t, y, P), 0:dt:T, y0, opts);
else
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    [t, y] = ode45(@(t, y) rhs(t, y, P), 0:dt:T, y0, opts);
end
nt = numel(t);
R = reshape(y(:, 1:5*N), nt, N, 5);
S.t = t;
S.rhoI = R(:,:,3) * f;
S.rhoR = R(:,:,4) * f;
end

%% ------------------------------------------------------------------------
function dy = rhs(~, y, P)
N = numel(P.f); f = P.f(:);
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

tot = sum(M, 2);
m   = sum(f .* tot);
Bi  = flipud(cumsum(flipud(f .* tot)));         % F_i b(t,h_i)
soc = (P.eta / P.tau) * ( P.p * (m * r - M) ...
      + (1 - P.p) * (bsxfun(@times, Bi, r) - bsxfun(@times, Fi, M)) );

dM = [ -(c + v) .* mS, ...
        c .* mS - P.sE * ME, ...
        P.sE * ME - P.gI * MI, ...
        P.gI * MI, ...
        v .* mS ] + soc;

dy = [dr(:); dM(:)];
end

%% ------------------------------------------------------------------------
function [I, R400] = montecarlo(P, seed, mode, T, dt, Na)
% Boltzmann model (B) with Na agents. States: 1=S, 2=E, 3=I, 4=R, 5=V.
% Opinions: rule w' = w + eta*chi*(w* - w), at rate 1/tau; chi = 1 if h* >= h
% and Bernoulli(p) if h* < h. Infection and vaccination at rates that depend on w.
rng(seed);
f = P.f(:)';
cf = cumsum(f);
u0 = rand(Na, 1);
lev = 1 + (u0 > cf(1)) + (u0 > cf(2));          % level of each agent (N = 3)
w = P.mS0(:); w = w(lev);
switch mode
    case 'unif'                                  % uniform spread inside each level
        dl = [0.2 1.0 0.1]; dl = dl(:);
        w = min(1, max(-1, w + dl(lev) .* (2*rand(Na,1) - 1)));
    case 'bimodal'                               % bimodal middle level at +-1
        mid = (lev == 2);
        w(mid) = 2*(rand(sum(mid),1) > 0.5) - 1;
end
st = ones(Na, 1);
st(rand(Na,1) < 0.02) = 3;                       % 2 % initially infectious
nst = round(T / dt);
I = zeros(nst + 1, 1);
I(1) = mean(st == 3);
for s = 1:nst
    % --- opinions ---
    idx = find(rand(Na,1) < dt / P.tau);
    j = randi(Na, numel(idx), 1);
    chi = (lev(j) >= lev(idx)) | (rand(numel(idx),1) < P.p);
    wn = w;
    wn(idx) = w(idx) + P.eta * chi .* (w(j) - w(idx));
    w = wn;
    % --- epidemic ---
    th = P.bT / 4 * sum(1 - w(st == 3)) / Na;
    pinf = 1 - exp(-dt * th * (1 - w));
    pvac = 1 - exp(-dt * P.bV / 2 * (1 + w));
    S_ = (st == 1);
    inf = S_ & (rand(Na,1) < pinf);
    vac = S_ & ~inf & (rand(Na,1) < pvac);
    new_st = st;
    new_st(inf) = 2;
    new_st(vac) = 5;
    new_st((st == 2) & (rand(Na,1) < 1 - exp(-dt * P.sE))) = 3;
    new_st((st == 3) & (rand(Na,1) < 1 - exp(-dt * P.gI))) = 4;
    st = new_st;
    I(s+1) = mean(st == 3);
end
R400 = mean(st == 4);
end
