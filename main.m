%% RP1
% Initalizations

% Constant
g = 9.8; % Gravity acelleration in m/s

% Equilibrium 
h = 1000;    % m
M = 0.20;
aa0 = 7.23 * pi /180; % rad
gg0 = 0;    % rad
u0 = 130.7 * 0.514444; % km/h
flaps = 20 * pi /180;  % rad
tt0 = aa0 + gg0;
w0 = 0; % ASSUMPTION

% Equilibrium Inputs
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
ybb = -0.0416; yp  = 0; yr  = 0; ydr = -0.013;
yda = 0; % Assumption
mwp = -0.0008;
lbb = -0.1211; lp = -0.5594; lr  = 0.0586; lda = -0.484; ldr = -0.024;
nbb = 0.2642; np  = -0.0028; nr  = -0.1134; nda = 0;ndr = -0.145;

% Computations
% Underscore represents the character '
l_bb = lbb + (Ixz/Ix) * nbb;
l_p  = lp  + (Ixz/Ix) * np;
l_r  = lr  + (Ixz/Ix) * nr; 
l_da = lda + (Ixz/Ix) * nda;
l_dr = ldr + (Ixz/Ix) * ndr;

n_bb = nbb + (Ixz/Iz) * lbb;
n_p  = np  + (Ixz/Iz) * lp;
n_r  = nr  + (Ixz/Iz) * lr;
n_da = nda + (Ixz/Iz) * lda;
n_dr = ndr + (Ixz/Iz) * ldr;

%% Defining the matrix A and B
A = zeros([4,4]);
% Line 1
A(1,1) = ybb; A(1,2) = (yp + w0)/u0;
A(1,3) = (yr - u0)/u0; A(1,4) = (g * cos(tt0))/u0;
% Line 2
A(2,1) = l_bb; A(2,2) = l_p; A(2,3) = l_r;
% Line 3
A(3,1) = n_bb; A(3,2) = n_p; A(3,3) = n_r;
% Line 4
A(4,2) = 1; A(4,3) = tan(tt0);

B = [yda ydr;
     l_da l_dr;
     n_da n_dr;
     0 0];

% Compute the poles, damping coefficients, frequencies, and time constants
damp(A)

%% Ponto 2
% Realimentação da razão de guinada para o rudder para establiziação do RH
rlocus(A, B(:,2),-[0 0 1 0], 0)
% k1 = rlocfind(A, B(:,2),-[0 0 1 0], 0)
k1 = 3.60;
A1 =  A + B(:,2)*[0 0 k1 0];
damp(A1)

% Falta estabilizar o rolamento puro
%% 
% Realimentação dos ailerons para a razão de rolamento para estabilizar o
% rolamento
rlocus(A1, B(:,1),-[0 1 0 0], 0)
% k2 = rlocfind(A1, B(:,1),-[0 1 0 0], 0)
k2 = 0.3377;
A_AF = A1 + B(:,1) * [0 k2 0 0];
damp(A_AF)




 
