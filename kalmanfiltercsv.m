%% Settings (tuned on FAR2026; re-tune for IREC)
data_file = 'FAR2026.csv';
t0 = 1230;  t1 = 1285;                 % analysis window, s
t_pad_end = 1240;                      % samples before this are "on the pad", s

T_pad = 32.75;                         % pad temperature, °C
h_pad = 628;                           % pad elevation, m above sea level

Rbaro  = 8.4;                          % baro noise variance, m^2
Raccel = 10;                           % accel noise variance (flight), (m/s^2)^2
Rgps   = 25;                           % GPS noise variance, m^2
sigma_j = 10;                          % process noise (jerk), m/s^3

v_trans_lo = 150;  v_trans_hi = 410;   % baro transonic gate, m/s
h_baro_max = 9144;                     % baro altitude gate (30k ft), m
v_gps_max  = 500;                      % GPS velocity limit, m/s
a_gps_max  = 39.24;                    % GPS acceleration limit (4 g), m/s^2

n_sigma     = 3;                       % innovation gate width
max_rej     = 25;                      % baro: accept after this many rejections
max_rej_gps = 5;                       % GPS: accept after this many rejections

%% Load data
data = readmatrix(data_file);
acc  = data(:, [1 2]);    acc  = acc(~any(isnan(acc), 2), :);
baro = data(:, [10 13]);  baro = baro(~any(isnan(baro), 2), :);
gps  = data(:, [48 51]);  gps  = gps(~any(isnan(gps), 2), :);

%% Prepare data
acc  = acc(acc(:,1)   >= t0 & acc(:,1)   <= t1, :);
baro = baro(baro(:,1) >= t0 & baro(:,1) <= t1, :);
gps  = gps(gps(:,1)   >= t0 & gps(:,1)   <= t1, :);
[~, iu] = unique(gps(:,1));  gps = gps(iu, :);       % drop repeated GPS rows

a_meas = (acc(:,2) - 1) * 9.81;                      % g -> m/s^2, remove 1 g

pad = baro(:,1) < t_pad_end;
baro(:,2) = baro(:,2) - mean(baro(pad, 2));          % zero baro to pad
baro_scale = (T_pad + 273.15) / (288.15 - 0.0065*h_pad);
baro(:,2) = baro(:,2) * baro_scale;                  % temperature correction

gpad = gps(:,1) < t_pad_end;
gps(:,2) = gps(:,2) - mean(gps(gpad, 2));            % zero GPS to pad

%% Filter setup
Hbaro  = [1 0 0];
Haccel = [0 0 1];
x_est = [0; 0; 0];
p_est = diag([1, 0.01, 0.1]);

Na = size(acc, 1);  Nb = size(baro, 1);  Ng = size(gps, 1);
j = 1;  g = 1;                                       % next baro / GPS sample
n_rej = 0;  n_rej_gps = 0;                           % consecutive rejections

est   = zeros(3, Na);  est(:,1) = x_est;
t_est = zeros(1, Na);  t_est(1) = acc(1,1);

%% Filter loop
for k = 2:Na
    % predict
    dt = acc(k,1) - acc(k-1,1);
    F  = [1 dt dt^2/2; 0 1 dt; 0 0 1];
    G  = [dt^3/6; dt^2/2; dt];
    Q  = sigma_j^2 * (G*G');
    x_pred = F * x_est;
    p_pred = F*p_est*F' + Q;

    % accel update
    K = (p_pred*Haccel')/(Haccel*p_pred*Haccel' + Raccel);
    x_pred = x_pred + K*(a_meas(k) - Haccel*x_pred);
    p_pred = (eye(3) - K*Haccel)*p_pred;

    % baro updates
    while j <= Nb && baro(j,1) <= acc(k,1)
        v = x_pred(2);
        if ~(v > v_trans_lo && v < v_trans_hi) && x_pred(1) < h_baro_max
            r = baro(j,2) - Hbaro*x_pred;
            S = Hbaro*p_pred*Hbaro' + Rbaro;
            if abs(r) <= n_sigma*sqrt(S) || n_rej >= max_rej
                K = (p_pred*Hbaro')/S;
                x_pred = x_pred + K*r;
                p_pred = (eye(3) - K*Hbaro)*p_pred;
                n_rej = 0;
            else
                n_rej = n_rej + 1;
            end
        end
        j = j + 1;
    end

    % GPS updates
    while g <= Ng && gps(g,1) <= acc(k,1)
        v = x_pred(2);  a = x_pred(3);
        if ~(abs(v) > v_gps_max || abs(a) > a_gps_max)
            r = gps(g,2) - Hbaro*x_pred;
            S = Hbaro*p_pred*Hbaro' + Rgps;
            if abs(r) <= n_sigma*sqrt(S) || n_rej_gps >= max_rej_gps
                K = (p_pred*Hbaro')/S;
                x_pred = x_pred + K*r;
                p_pred = (eye(3) - K*Hbaro)*p_pred;
                n_rej_gps = 0;
            else
                n_rej_gps = n_rej_gps + 1;
            end
        end
        g = g + 1;
    end

    % finish step
    x_est = x_pred;
    p_est = p_pred;
    est(:,k) = x_est;
    t_est(k) = acc(k,1);
end

%% Results
[h_max, i_h] = max(est(1,:));  fprintf('apogee %.0f m at %.2f s\n', h_max, t_est(i_h));
[v_max, i_v] = max(est(2,:));  fprintf('max velocity %.0f m/s at %.2f s\n', v_max, t_est(i_v));

figure; plot(t_est, est(1,:)); hold on;
plot(baro(:,1), baro(:,2)); plot(gps(:,1), gps(:,2));
legend('estimate', 'baro', 'GPS'); title('Altitude'); xlabel('time (s)'); xlim([t0 t1]);

figure; plot(t_est, est(2,:)); title('Velocity'); xlabel('time (s)'); xlim([t0 t1]);

figure; plot(t_est, est(3,:)); hold on; plot(acc(:,1), a_meas);
legend('estimate', 'accel'); title('Acceleration'); xlabel('time (s)'); xlim([t0 t1]);