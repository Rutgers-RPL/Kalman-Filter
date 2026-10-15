x_true = 100;
N = 1000;

x_est = 0;
p_est = 10000;
q = .1;
r = 25;

est = zeros(1, N);
meas = zeros(1, N);
gains = zeros(1,N);

for i = 1:N
    z = x_true + 5*randn;            
    x_pred = x_est;
    p_pred = p_est + q;
    K = (p_pred)/(p_pred+r);
    x_est = x_pred + K*(z - x_pred); 
    p_est = (1-K)*p_pred;
    est(i) = x_est;
    meas(i) = z;
    gains(i) = K;
end

plot(meas)
hold on
plot(est)
yline(x_true)
legend('measured', 'estimate', 'truth')
figure;
plot(gains);