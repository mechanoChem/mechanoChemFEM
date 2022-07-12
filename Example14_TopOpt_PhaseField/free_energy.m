c = 0.0:0.0001:1.0;
A=1.0; B=3.5*A;
phi=A*(c.*log(c) + (1-c).*log(1-c)) + B*c.*(1-c);
mu = A*log(c./(1-c))+B*(1-2*c);

plot(c,phi)
plot(c,mu)
