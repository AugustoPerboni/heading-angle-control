clear all
close all
clc

%% RP1
% Inicializaçoes
states_seguimento = {'\beta' 'p' 'r' '\phi'};
inputs_seguimento = {'delta_a' 'delta_r'};
outputs_seguimento  = {'vel_horiz' 'vel_vert' 'pitch rate' 'vel_subida' 'altitude'};

% Constantes
g = 9.8; % Gravity acelleration in m/s

% Equilibrio
h = 1000;    % m
M = 0.20;
aa0 = 7.23 * pi /180; % rad
gg0 = 0;    % rad
u0 = 130.7 /1.94384; % km/h
flaps = 20 * pi /180;  % rad
tt0 = aa0 + gg0;
w0 = 0; % ASSUMPTION

% Equilibrio Inputs
da0 = 0.00; % rad
de0 = 0.00; % rad
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
xu=-0.0745; xw=0.0002; 
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
A = [ybb   yp+w0/u0 yr - u0/u0 (g * cos(tt0))/u0; ...
            l_bb   l_p       l_r           0; ...
            n_bb   n_p       n_r           0; ...
            0       1        tan(tt0)       0]

B = [yda ydr;
     l_da l_dr;
     n_da n_dr;
     0 0];
C = eye(4);
D = zeros(4,2);

% Compute the poles, damping coefficients, frequencies, and time constants
damp(A)

%% Ponto 2
% Realimentação da razão de guinada 
sys_1 = ss(A, B(:,2),[0 0 1 0], 0);
figure;
rlocus(-sys_1)
k1 = 5.25;
A1 =  A - B(:,2)*[0  0 -k1 0];
damp(A1)
figure;
step(ss(A1,B,C,D,'inputname', inputs_seguimento, 'outputname' , {'\beta' 'p' 'r' '\phi'}))
title('Resposta do sistema');
xlabel('Tempo');
ylabel('Amplitude');

% Falta estabilizar o rolamento puro
%% 
% Realimentação dos do rolamento
sys_2=ss(A1, B(:,1),[0 1 0 0], 0);
figure;
rlocus(-sys_2)
% rlocfind(-sys_2)

% k2 = rlocfind(A1, B(:,1),-[0 1 0 0], 0)
k2 = 0.41;
A_AF = A1 + B(:,1) * [0 k2 0 0];
damp(A_AF)
C = eye(4);
D = zeros(4,2);
sys_3=ss(A_AF, B,C, D);
figure;
step(ss(A_AF,B,C,D,'inputname', inputs_seguimento, 'outputname' , {'\beta' 'p' 'r' '\phi'}))
title('Resposta do sistema');
xlabel('Tempo');
ylabel('Amplitude');
%%

t=0:0.1:100;
% Degray de 15 graus
u=zeros(1001,1);
u(30:1001 , 1) = deg2rad(15);
% u=(-heaviside(t-15)+1)*deg2rad(15); 
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
A_int = [A zeros(4,1); zeros(1,5)];
A_int(5,3) = 1/cos(tt0);
B_int = [B; 0 0];
C_int = eye(5);
D_int = zeros(5,2);

% Valores maximos
damax = deg2rad(30);        % rad
drmax = deg2rad(30);        % rad

bbmax = deg2rad(15);        % rad
phimax = deg2rad(30);        % rad ( angulo de rolamento )
pmax = 0.025;
rmax = 0.01;
psimax = deg2rad(60);


% Sem estados integrativos
% R = 0.3*diag([1/(damax*0.1)^2 2*1/(drmax*0.1)^2]); % 10% do max 
% Q = diag([5*1/bbmax^2 0.23*1/pmax^2 1/rmax^2 5*1/phimax^2 20*1/psimax^2]);
% Klqr = lqr(A_int,B_int,Q,R);
% A_AF_int = A_int - B_int*Klqr;
% sys_lqr = ss(A_AF_int, B_int, C_int, D_int);
% k_sys_lqr = dcgain(sys_lqr);
% prek_sys_lqr = pinv([k_sys_lqr(1,1), k_sys_lqr(1,2);k_sys_lqr(5,1), k_sys_lqr(5,2)]);
% damp(A_AF);


% Adicionando estados integrativos
Ai = [A_int zeros(5,2); 1 0 0 0 0 0 0; 0 0 0 0 1 0 0];
Bi = [B_int; 0 0; 0 0];
Ci = eye(7);
Di = zeros(7,2);

% xm2= [bbmax pmax rmax phimax lambdamax];
xm2= [1 4 4 10 50 1 1];
um2 = [damax drmax];       
Qi = diag(xm2);%.^(-2)
Ri = diag(um2.^(-2));
Klqri = lqr(Ai, Bi, Qi, Ri);
Ai_AF = Ai - Bi * Klqri;
damp(Ai_AF)






% % Valores maximos
% damax = deg2rad(30);        % rad
% drmax = deg2rad(30);        % rad
% bbmax = deg2rad(3);        % rad
% phimax = deg2rad(30);        % rad
% pmax = deg2rad(20);
% rmax = deg2rad(20);
% psimax = deg2rad(40);
% lambdamax= deg2rad(180);
% 
% % xm2= [bbmax pmax rmax phimax lambdamax];
% xm2= [1 4 4 10 20];
% um2 = [damax drmax];       
% Q2 = diag(xm2);  %.^(-2)
% R2 = diag(um2.^(-2));
% 
% % adição dos estado integrador para controlo da atitude
% A_int = [A zeros(4,1); zeros(1,5)];
% A_int(5,4) = g/u0;
% A_int
% damp(A_int)
% B_int = [B; 0 0;];
% B_int
% C_int = eye(5);
% D_int = zeros(5,2);
% 
% 
% 
% klqr_int = lqr(A_int,B_int,Q2,R2)
% AF_int = A_int - B_int*klqr_int;
% damp(AF_int);
% sys_lqr=ss(AF_int, B_int, C_int, D_int);
% k_sys_lqr=dcgain(sys_lqr)
% prek_sys_lqr=pinv([k_sys_lqr(5,1) k_sys_lqr(5,2)])
% 
% out=sim("Ponto3.slx")
% 
% figure;
% hold on
% grid on
% plot(out.lambda)
% plot(out.p)
% plot(out.r)
% plot(out.phi)
% plot(out.beta)
% title('Resposta das variáveis de estado')
% legend('\lambda','p','r','\phi','\beta')
% xlabel('tempo(segundos)')
% ylabel('amplitude(rad ou rad/s)')


