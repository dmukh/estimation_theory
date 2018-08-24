function Xdot = populateXdot(t, X)

global re h0 gamma0 rho0 omegaE 

phi_i = reshape(X(13:end),12,12);
X_i = X(1:12);
x = X_i(1); y = X_i(2); z = X_i(3); xdot = X_i(4); ydot = X_i(5); zdot = X_i(6);
mu = X_i(7); J2 = X_i(8); beta = X_i(9); xs = X_i(10); ys = X_i(11); zs = X_i(12);

a = xdot + y*omegaE;
b = ydot - x*omegaE;
c = zdot;
Q = norm([a b c]);
r = sqrt(x^2 + y^2 + z^2);
rho = rho0*exp(-gamma0*(r-re-h0));

L = -mu*x/r^3*(1 - 3/2*(re/r)^2*J2*(5*(z/r)^2-1));
M = -mu*y/r^3*(1 - 3/2*(re/r)^2*J2*(5*(z/r)^2-1));
N = -mu*z/r^3*(1 - 3/2*(re/r)^2*J2*(5*(z/r)^2-3));


Xdot = zeros(156,1);

Xdot(1:12) = [xdot; ydot; zdot; L-beta*rho*Q*a; M-beta*rho*Q*b; N-beta*rho*Q*c; 0;0;0;0;0;0];

A = populateA(X_i);
phidot = A*phi_i;
Xdot(13:end) = reshape(phidot, 144, 1);





