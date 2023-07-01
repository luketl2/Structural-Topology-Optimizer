Fax = 0; % N
Fsh = 1000; % N
L = 200; % mm
M = Fsh*L; % N*mm
t = 10; % mm
h = 100; % mm
cMax = h/2; % mm
A = h*t; % mm^2
sigmaAxial = Fax/A; % MPa
tao = Fsh/A; % MPa
Ix = t*h^3/12; % mm^4
sigmaBending = M*cMax/Ix; % MPa
sigmaVM = sqrt((sigmaAxial+sigmaBending)^2+3*tao^2) % MPa