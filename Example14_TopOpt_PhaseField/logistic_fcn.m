clear all
close all

x0=0.5;
A=1;
B=10;

x=-0.1:0.01:1.1;
y=A./(1.0+exp(-B .* (x-x0)));

dydx= 1484.13*exp(10*x) ./(148.413+exp(10*x)).^2;

figure()
plot(x, y)
figure()
plot(x, dydx)