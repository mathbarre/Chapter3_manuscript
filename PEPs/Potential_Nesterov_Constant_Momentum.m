
% Code for solving the feasibility problem (3.35)
clc; clear all;

verbose   = 1;
tolerance = 1e-8;


L  = 2;
mu = 0.1;

beta = (sqrt(L)-sqrt(mu))/(sqrt(L)+sqrt(mu));

dimG = 5; dimF = 3;
x0 = zeros(dimG,1);x0(1)=1;
x1 = zeros(dimG,1);x1(2)=1;


g0 = zeros(dimG,1);g0(3)=1;
g1 = zeros(dimG,1);g1(4)=1;
g2 = zeros(dimG,1);g2(5)=1;

x2 = (1+beta)*x1 - (1+beta)/L*g1 -beta*x0 + beta/L*g0;

gs = zeros(dimG,1);
xs = zeros(dimG,1);


f0  = zeros(dimF,1); f0(1) = 1;
f1  = zeros(dimF,1); f1(2) = 1;
f2  = zeros(dimF,1); f2(3) = 1;
fs  = zeros(dimF,1);

nbPts = 4;
XPEP  = [xs x0 x1 x2];
GPEP  = [gs g0 g1 g2];
FPEP  = [fs f0 f1 f2];

% MATRICES




rho = (1-sqrt(mu/L));


S = sdpvar(3);
% the potential we found
%S = L^2/(2*(L-mu))*[(1+sqrt(mu/L))^2, - (1+sqrt(mu/L)), (1+sqrt(mu/L));-(1+sqrt(mu/L)),1,-1;(1+sqrt(mu/L)),-1,1];

cons_SDP = [x2-xs x1-xs g1/L]*S*[x2-xs x1-xs g1/L].' - rho*[x1-xs x0-xs g0/L]*S*[x1-xs x0-xs g0/L].';


cons_LIN = (f1-fs) - rho*(f0-fs);

nu = sdpvar(nbPts,nbPts,'full');

for i = 1:nbPts
    for j = 1:nbPts
        if j ~= i %&& (i~=4) && (j~=4) && i~=1 && i~=2
            xi = XPEP(:,i); xj = XPEP(:,j);
            gi = GPEP(:,i); gj = GPEP(:,j);
            fi = FPEP(:,i); fj = FPEP(:,j);

            
            cons_SDP = cons_SDP - nu(i,j) * (gi*(xj-xi).' + 1/2/(L-mu)*(gi-gj- mu*(xi-xj))*(gi-gj - mu*(xi-xj)).' ...
                       + mu/2*(xi-xj)*(xi-xj).');
            cons_LIN = cons_LIN - nu(i,j) * (fi - fj);
        end
    end
end




size(cons_SDP)
cons_SDP = (cons_SDP.' + cons_SDP)/2;
cons = cons_SDP <=0;
cons = cons+ (cons_LIN == 0);
cons = cons + (S >= 0);
cons = cons + (nu >= 0);



% in order to have simple sum of squares residual
obj =trace(-cons_SDP);



solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
solverDetails=optimize(cons,obj,solver_opt);