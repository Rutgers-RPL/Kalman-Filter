%data import
data = readmatrix('FAR2026.csv');

acc = data(:, [1 2]);
acc = acc(~any(isnan(acc), 2), :);

baro = data(:, [10 13]);
baro = baro(~any(isnan(baro), 2), :);

gps = data(:, [48 51]);
gps = gps(~any(isnan(gps), 2), :);




%initial flight data
t0 = 1230;  t1 = 1285;
acc  = acc(acc(:,1)   >= t0 & acc(:,1)   <= t1, :);
baro = baro(baro(:,1) >= t0 & baro(:,1) <= t1, :);
gps  = gps(gps(:,1)   >= t0 & gps(:,1)   <= t1, :);

a_meas = (acc(:,2) - 1) * 9.81;

pad = baro(:,1) < 1240;
baro(:,2) = baro(:,2) - mean(baro(pad, 2));
gpad = gps(:,1) < 1240;
gps(:,2) = gps(:,2) - mean(gps(gpad, 2));





Rbaro   = 8.4;
Raccel  = 10;
sigma_j = 10;
Hbaro   = [1 0 0];
Haccel  = [0 0 1];

x_est = [0; 0; 0];
p_est = diag([1, 0.01, 0.1]);

Na = size(acc, 1);      % number of accel samples
Nb = size(baro, 1);     % number of baro samples
j  = 1;                 % pointer: next baro sample to use

est   = zeros(3, Na);   % one column of [h; v; a] per step
t_est = zeros(1, Na);



est(:,1) = x_est;
t_est(1) = acc(1,1);
n_rej = 0;          % consecutive rejections
max_rej = 25;       % then accept anyway
% Initialize the first estimate and timestamp from the initial measurements
for k = 2:Na
    dt = acc(k,1) - acc(k-1,1);
    F  = [1 dt dt^2/2; 0 1 dt; 0 0 1];
    G  = [dt^3/6; dt^2/2; dt];
    Q  = sigma_j^2 * (G*G');

    x_pred = F * x_est;
    p_pred = F*p_est*F' + Q;

    K = (p_pred*Haccel')/(Haccel*p_pred*Haccel' + Raccel);
    x_pred = x_pred + K*(a_meas(k) - Haccel*x_pred);
    p_pred = (eye(3) - K*Haccel)*p_pred;

    while j <= Nb && baro(j,1) <= acc(k,1)
        v = x_pred(2);
        v = x_pred(2);
        if ~(v > 150 && v < 410)                       % velocity gate 
            r = baro(j,2) - Hbaro*x_pred;              % innovation
            S = Hbaro*p_pred*Hbaro' + Rbaro;           % expected size
            if abs(r) <= 3*sqrt(S) || n_rej >= max_rej % innovation gate, with escape
                K = (p_pred*Hbaro')/S;
                x_pred = x_pred + K*r;
                p_pred = (eye(3) - K*Hbaro)*p_pred;
                n_rej = 0;                             % reset the counter
            else
                n_rej = n_rej + 1;                     % count the rejection
            end
        end
        j = j + 1;
        
    end

    x_est = x_pred;
    p_est = p_pred;
    est(:,k) = x_est;
    t_est(k) = acc(k,1)
end

[h_max, i_h] = max(est(1,:));   fprintf('apogee %.0f m at %.2f s\n', h_max, t_est(i_h));
[v_max, i_v] = max(est(2,:));   fprintf('max velocity %.0f m/s at %.2f s\n', v_max, t_est(i_v));

%plots and data blah blah blah
figure; plot(t_est, est(1,:)); hold on;
plot(baro(:,1), baro(:,2)); plot(gps(:,1), gps(:,2));
legend('estimate', 'baro', 'GPS'); title('Altitude'); xlabel('time (s)');

figure; plot(t_est, est(2,:)); title('Velocity'); xlabel('time (s)');

figure; plot(t_est, est(3,:)); hold on; plot(acc(:,1), a_meas);
legend('estimate', 'accel'); title('Acceleration'); xlabel('time (s)');

xlim([1240 1285])