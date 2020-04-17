function i= int_torque(a, r)
a0= a(1:end-1);
a1= a(2:end);
dR= diff(r);
R0= r(1:end-1);

i= ((2*a1+a0).*dR.^2 + (3*R0.*a1 + 3*R0.*a0).*dR)/6;
