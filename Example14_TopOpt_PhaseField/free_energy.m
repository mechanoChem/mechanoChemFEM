close all
clear all
c = 0.0:0.0001:1.0;

beta1 = 0.5;
beta2 = 50;
phi=c.^2 .* (1-c).^2 + beta1 * (power(10, -beta2 * c) + power(10, beta2 * (c-1)));

% already hard code beta1 and beta2 into the chemical potential
mu = 4 * c.^3 - 6*c.^2 + 2 .*c - 57.5646 * power(10, -50*c) + 5.75646*power(10, -49) * power(10, 50*c);

figure()
plot(c,phi)
figure()
plot(c,mu)
