%% Inicializations

% Constant
g = 9.8; % Gravity acelleration in m/s

% Equilibrium 
h = 1000;    % m
M = 0.20;
aa0 = 7.23 * pi /180; % rad
gg0 = 8 * pi /180;    % rad
u0 = 130.7 ; % km/h
flaps = 20 * pi /180;  % rad
theta0 = 0; % ASSUMPTION
w0 = 0; % ASSUMPTION

% Equilibrium Inputs
th0 = 0.57; % percentage
de0 = 0.00; % rad
da0 = 0.00; % rad
dr0 = 0.00; % rad 

% Maximum values
Teng = 5.00;       % s
demax = +26/ -26 * pi /180 ; % rad
damax = 30 * pi /180;        % rad
drmax = 30 * pi /180;        % rad
flapmax = 40 * pi /180;      % rad

% Inertial data
m = 49051;    % kg 
Ix = 3095873; % kg . m ^2
Iy = 4303047; % kg . m ^2
Iz = 4526287; % kg . m ^2
Ixz = 1356;   % kg . m ^2

% Wing data 
S = 105.63; % m ^2
b = 28.880; % m 
c = 3.658;  % m
aamax = 13.12 * pi /180; % rad

% Derivatives ( no units or SI units)

xu  = -0.0745;
xw  = 0.002;
xde = 0;
xdf = -0.606;
xdt = 4.088;

ybb = -0.0416;
yp  = 0;
yr  = 0;
Ydr = -0.013;

zu  = -0.2917;
zw  = -0.6220;
zwp = -0.0139;
zq  = -2.5046;
zde = -1.349;
zdf = -8.128;
zdt = 0;

mu  = 0;
mw  = -0.0146;
mq  = -0.1230;
mwp = -0.0008;
mde = -0.459;
mdf = 0.002;
mdt = 0.036;


lbb = -0.1211; 
lp  = -0.5594;
lr  = 0.0586;
Lda = -0.484;
Ldr = -0.024;

nbb = 0.2642;
np  = -0.0028;
nr  = -0.1134;
Nda = 0;
Ndr = -0.145;

% Computations
% Underscore represents the character '

% l_v = lv + (Ixz/Ix)*nv 
l_p = lp + (Ixz/Ix) * np;
l_r = lr + (Ixz/Ix) * nr; 

% n_v = nv + (Ixz/Iz) * lv;
n_p = np + (Ixz/Iz) * lp;
n_r = nr + (Ixz/Iz) * lr;


%% Defining the matrix A

A = zeros([5,5]);

% Line 1
% A(1,1) = yv
A(1,2) = yp + w0;
A(1,3) = yr - u0;
A(1,4) = g * cos(theta0);

% Line 2
% A(2,1) = l_v;
A(2,2) = l_p;
A(2,3) = l_r;

% Line 3
% A(3,1) = l_v;
A(3,2) = n_p;
A(3,3) = n_r;

% Line 4
A(4,2) = 1;
A(4,3) = tan(theta0);

% Line 5 
A(5,3) = 1/cos(theta0);



 
