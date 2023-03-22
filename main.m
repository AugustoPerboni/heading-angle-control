%% RP1
% Inicializaçoes

% Constantes
g = 9.8; % Gravity acelleration in m/s

% Equilibrio
h = 1000;    % m
M = 0.20;
aa0 = 7.23 * pi /180; % rad
gg0 = 0;    % rad
u0 = 130.7 * 0.514444; % km/h
flaps = 20 * pi /180;  % rad
tt0 = aa0 + gg0;
w0 = 0; % ASSUMPTION

% Equilibrio Inputs
da0 = 0.00; % rad
dr0 = 0.00; % rad 

% Dados inerciais
m = 49051;    % kg 
Ix = 3095873; % kg . m ^2
Iy = 4303047; % kg . m ^2
Iz = 4526287; % kg . m ^2
Ixz = 1356;   % kg . m ^2

% Dados da asa
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

% Defining the matrix A and B
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
% Realimentação da razão de guinada 
rlocus(A, B(:,2),-[0 0 1 0], 0)
% k1 = rlocfind(A, B(:,2),-[0 0 1 0], 0)
k1 = 4.15;
A1 =  A + B(:,2)*[0 0 k1 0];
damp(A1)

% Falta estabilizar o rolamento puro
%% 
% Realimentação dos do rolamento
rlocus(A1, B(:,1),-[0 1 0 0], 0)
% k2 = rlocfind(A1, B(:,1),-[0 1 0 0], 0)
k2 = 0.3377;
A_AF = A1 + B(:,1) * [0 k2 0 0];
damp(A_AF)
%%

t=0:0.001:100;
% Degray de 15 graus
u=(-heaviside(t-15)+1)*deg2rad(15); 
plot(u)

% Definição dos sistemas a serem plotados
C = eye(4);
D = zeros(4,2);
sys_aa=ss(A,B(:,1),eye(4),D(:,1),'OutputName',{'\beta','p','r','\phi'}); % ailerons (sem realimentação)
sys_ra=ss(A,B(:,2),eye(4),D(:,2),'OutputName',{'\beta','p','r','\phi'}); % rudder (...)
sys_af=ss(A_AF,B(:,1),eye(4),D(:,1),'OutputName',{'\beta','p','r','\phi'}); % ailerons (com realimentação)
sys_rf=ss(A_AF,B(:,2),eye(4),D(:,2),'OutputName',{'\beta','p','r','\phi'}); % rudder (...)

% Sistema nao estabilizado
figure();
lsim(sys_aa,u,t); 
title('Anel aberto - Degrau ailerons')
figure();
lsim(sys_ra,u,t); 
title('Anel aberto - Degrau rudder')
% Sistema estabilizado
figure();
lsim(sys_af,u,t); 
title('Sistema estabilizado - Degrau ailerons')
figure();
lsim(sys_rf,u,t); 
title('Sistema estabilizado - Degrau rudder')

%% Ponto 3
% A matriz A passa a ter o angulo de guinada como variavel de estado
A = [A zeros(4,1); zeros(1,5)];
A(5,3) = 1/cos(tt0);
B = [B; 0 0];
C = eye(5);
D = zeros(5,2);

% Valores maximos
damax = deg2rad(30);        % rad
drmax = deg2rad(30);        % rad

bbmax = deg2rad(15);        % rad
phimax = deg2rad(30);        % rad ( angulo de rolamento )
pmax = 0.025;
rmax = 0.01;
psimax = deg2rad(60);

R = 0.3*diag([1/(damax*0.1)^2 2*1/(drmax*0.1)^2]); % 10% do max 
Q = diag([5*1/bbmax^2 0.23*1/pmax^2 1/rmax^2 5*1/phimax^2 20*1/psimax^2]);

Klqr = lqr(A,B,Q,R);

A_AF = A - B*Klqr;
K = dcgain(A_AF,B,C,D);
K_pre = [K(1,1), K(1,2);K(5,1), K(5,2)];
damp(A_AF);


close all;











 
