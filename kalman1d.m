N = 1000;

dt = .1;
x_est = [0; 0];
p_est = 10000 * eye(2);
Q = [.01 0; 0 1];
R = 25;
F = [1 dt; 0 1];
H = [1 0];
truth = zeros(1, N);
truth_v = zeros(1, N);
est_h = zeros(1, N);
est_v = zeros(1, N);
meas = zeros(1, N);
gains = zeros(1, N);
for i = 1:N
    x_true = .5*20*(i*dt)^2;
    z = x_true + 5*randn;            
    x_pred = F * x_est;
    p_pred = F*p_est*F' + Q;
    K = (p_pred*H')/(H*p_pred*H'+R);
    x_est = x_pred + K*(z - H* x_pred); 
    p_est = (eye(2)-K * H)*p_pred;
    truth(i) = x_true;
    est_h(i) = x_est(1);
    est_v(i) = x_est(2);
    meas(i) = z;
    gains(i) = K(1);
    truth_v(i) = 20 * i * dt;
end

plot(meas); hold on; plot(est_h); plot(truth);
legend('measured', 'estimate', 'truth'); title('Altitude');

figure; plot(est_v); hold on; plot(truth_v);
legend('estimate', 'truth'); title('Velocity');

figure; plot(gains); title('K altitude');
figure; plot(est_h - truth); title('Altitude error');