function A = populateA(X)

global re h0 gamma0 rho0 omegaE 

x = X(1); y = X(2); z = X(3); xdot = X(4); ydot = X(5); zdot = X(6);
mu = X(7); J2 = X(8); beta = X(9); xs = X(10); ys = X(11); zs = X(12);

a = xdot + y*omegaE;
b = ydot - x*omegaE;
c = zdot;
Q = norm([a b c]);
r = sqrt(x^2 + y^2 + z^2);
rho = rho0*exp(-gamma0*(r-re-h0));

L = -mu*x/r^3*(1 - 3/2*(re/r)^2*J2*(5*(z/r)^2-1));
M = -mu*y/r^3*(1 - 3/2*(re/r)^2*J2*(5*(z/r)^2-1));
N = -mu*z/r^3*(1 - 3/2*(re/r)^2*J2*(5*(z/r)^2-3));


f4X1 = L/x + mu*x^2/r^5*(3-15/2*(re/r)^2*J2*(7*z^2/r^2-1)) + beta*rho*a*(gamma0*Q*x/r + b*omegaE/Q);
f4X2 = mu*x*y/r^5*(3-15/2*(re/r)^2*J2*(7*z^2/r^2-1)) - beta*rho*(-gamma0*Q*a*y/r + a^2*omegaE/Q + Q*omegaE);
f4X3 = mu*x*z/r^5*(3-15/2*(re/r)^2*J2*(7*z^2/r^2-3)) + beta*rho*a*gamma0*Q*z/r;
f4X4 = -beta*rho*(a^2/Q+Q);
f4X5 = -beta*rho*a*b/Q;
f4X6 = -beta*rho*a*c/Q;
f4X7 = L/mu;
f4X8 = 3/2*(re/r)^2*(mu*x/r^3)*(5*z^2/r^2-1);
f4X9 = -rho*Q*a;

f5X1 = mu*y*x/r^5*(3-15/2*(re/r)^2*J2*(7*z^2/r^2-1)) + beta*rho*(gamma0*Q*b*x/r + b*omegaE/Q + Q*omegaE);
f5X2 = M/y + mu*y^2/r^5*(3-15/2*(re/r)^2*J2*(7*z^2/r^2-1)) - beta*rho*b*(-gamma0*Q*y/r + a*omegaE/Q);
f5X3 = mu*y*z/r^5*(3-15/2*(re/r)^2*J2*(7*z^2/r^2-3)) + beta*rho*b*Q*gamma0*z/r;
f5X4 = -beta*rho*a*b/Q;
f5X5 = -beta*rho*(b^2/Q + Q);
f5X6 = -beta*rho*b*c/Q;
f5X7 = M/mu;
f5X8 = 3/2*(re/r)^2*(mu*y/r^3)*(5*z^2/r^2-1);
f5X9 = -rho*Q*b;

f6X1 = mu*z*x/r^5*(3-15/2*(re/r)^2*J2*(7*z^2/r^2-3)) + beta*rho*c*(gamma0*Q*x/r + b*omegaE/Q);
f6X2 = mu*z*y/r^5*(3-15/2*(re/r)^2*J2*(7*z^2/r^2-3)) - beta*rho*c*(-gamma0*Q*y/r + a*omegaE/Q);
f6X3 = N/z + mu*z^2/r^5*(3-15/2*(re/r)^2*J2*(7*z^2/r^2-5)) + beta*rho*c*Q*gamma0*z/r;
f6X4 = -beta*rho*a*c/Q;
f6X5 = -beta*rho*b*c/Q;
f6X6 = -beta*rho*(c^2/Q + Q);
f6X7 = N/mu;
f6X8 = 3/2*(re/r)^2*(mu*z/r^3)*(5*z^2/r^2-3);
f6X9 = -rho*Q*c;

% A = dFdX
A = [0 0 0 1 0 0 0 0 0 0 0 0;...
    0 0 0 0 1 0 0 0 0 0 0 0;...
    0 0 0 0 0 1 0 0 0 0 0 0;...
    f4X1 f4X2 f4X3 f4X4 f4X5 f4X6 f4X7 f4X8 f4X9 0 0 0;...
    f5X1 f5X2 f5X3 f5X4 f5X5 f5X6 f5X7 f5X8 f5X9 0 0 0;...
    f6X1 f6X2 f6X3 f6X4 f6X5 f6X6 f6X7 f6X8 f6X9 0 0 0;...
    zeros(6, 12)];



