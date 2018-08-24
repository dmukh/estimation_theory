function [Htilde, G] = populateMeas(t, X)

r = X(1:3);

% Station locations
rs(1,:) = X(10:12);
rs(2,:) = [-2428826.1117 -4799750.4339 3417273.0738];
rs(3,:) = [-1736003.0850 -4425049.6149 4241427.1084];

global omegaE theta0

theta = theta0 + omegaE*t;

R = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];

rsI(1,:) = R*rs(1,:)';
rsI(2,:) = R*rs(2,:)';
rsI(3,:) = R*rs(3,:)';

% Populate G
Rsi = [r;r;r] - rsI;
G = [norm(Rsi(1,:)); norm(Rsi(2,:)); norm(Rsi(3,:))];

% Populate Htilde
Htilde = zeros(3,12);

Htilde(1,1) = Rsi(1,1)/norm(Rsi(1,:));
Htilde(1,2) = Rsi(1,2)/norm(Rsi(1,:));
Htilde(1,3) = Rsi(1,3)/norm(Rsi(1,:));
Htilde(1,10) = -(1/norm(Rsi(1,:)))*(Rsi(1,1)*cos(theta) + Rsi(1,2)*sin(theta));
Htilde(1,11) = -(1/norm(Rsi(1,:)))*(-Rsi(1,1)*sin(theta) + Rsi(1,2)*cos(theta));
Htilde(1,12) = -(1/norm(Rsi(1,:)))*Rsi(1,3);

Htilde(2,1) = Rsi(2,1)/norm(Rsi(2,:));
Htilde(2,2) = Rsi(2,2)/norm(Rsi(2,:));
Htilde(2,3) = Rsi(2,3)/norm(Rsi(2,:));

Htilde(3,1) = Rsi(3,1)/norm(Rsi(3,:));
Htilde(3,2) = Rsi(3,2)/norm(Rsi(3,:));
Htilde(3,3) = Rsi(3,3)/norm(Rsi(3,:));


