N = 150;
h_t= 0;
v_t= 0;
dt = .1;
x_est = [0; 0; 0];
p_est = 10000 * eye(3);
G = [dt^3/6; dt^2/2; dt];
sigma_j = 10;
Q = sigma_j^2 * (G*G');
Rbaro = 8.4;
Raccel = .006;
F = [1 dt (1/2)*dt^2; 0 1 dt;0 0 1];
Hbaro = [1 0 0];
Haccel = [0 0 1];
truth = zeros(1, N);
truth_v = zeros(1, N);
truth_a = zeros(1,N);
est_h = zeros(1, N);z
est_v = zeros(1, N);
est_a = zeros(1, N);
meas = zeros(1, N);
gains = zeros(1, N);
bound = zeros(1, N);
for i = 1:N
    t = i*dt;
    if t < 3
        a_t = 20;
    else 
        a_t = -9.81;
    end
    h_t = h_t + v_t*dt + 0.5*a_t*dt^2;
    v_t = v_t + a_t*dt;
    z_baro  = h_t + sqrt(Rbaro)*randn;
    z_accel = a_t + sqrt(Raccel)*randn;        
    x_pred = F * x_est;
    p_pred = F*p_est*F' + Q;
    K = (p_pred*Haccel')/(Haccel*p_pred*Haccel'+Raccel);
    x_pred = x_pred + K*(z_accel - Haccel* x_pred); 
    p_pred = (eye(3)-K * Haccel)*p_pred;
    K = (p_pred*Hbaro')/(Hbaro*p_pred*Hbaro' + Rbaro);
    x_est = x_pred + K*(z_baro - Hbaro*x_pred);
    p_est = (eye(3) - K*Hbaro)*p_pred;
    truth(i) = h_t;
    est_h(i) = x_est(1);
    est_v(i) = x_est(2);
    meas(i) = z_baro;
    gains(i) = K(1);
    truth_v(i) = v_t;
    est_a(i) = x_est(3);
    truth_a(i) = a_t;
    bound(i) = 3*sqrt(p_est(1,1));
end

plot(meas); hold on; plot(est_h); plot(truth);
legend('measured', 'estimate', 'truth'); title('Altitude');

figure; plot(est_v); hold on; plot(truth_v);
legend('estimate', 'truth'); title('Velocity');

figure; plot(gains); title('K altitude');
figure; plot(est_h - truth); title('Altitude error');

figure; plot(est_a); hold on; plot(truth_a);
legend('estimate', 'truth'); title('Acceleration');

figure; plot(est_h - truth); hold on; plot(bound); plot(-bound);
legend('error', '+3\sigma', '-3\sigma'); title('Altitude error vs bounds');