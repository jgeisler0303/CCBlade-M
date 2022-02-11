function i= int_torque_m(a, r)
a0= a(:, 1:end-1);
a1= a(:, 2:end);
dR= repmat(diff(r)', size(a, 1), 1);
R0= repmat(r(1:end-1)', size(a, 1), 1);

i= ((2*a1+a0).*dR.^2 + (3*R0.*a1 + 3*R0.*a0).*dR)/6;
