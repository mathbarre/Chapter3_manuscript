% Code for solving the feasibility problem (3.41)
clc; clear all;

verbose   = 1;
tolerance = 1e-6;


L=1;

alpha= 1.2;
lambda = alpha/L;


dimG = 3; dimF = 3;
x0 = zeros(dimG,1);x0(1)=1;

g0 = zeros(dimG,1);g0(2)=1;
g1 = zeros(dimG,1);g1(3)=1;

% gradient update
x1 = x0 - lambda*g0;

gs = zeros(dimG,1);
xs = zeros(dimG,1);


f0  = zeros(dimF,1); f0(1) = 1;
f1  = zeros(dimF,1); f1(2) = 1;

fs  = zeros(dimF,1);

nbPts = 3;
XPEP  = [xs x0 x1];
GPEP  = [gs g0 g1];
FPEP  = [fs f0 f1];



% set a value of t0
t0 = 2;
t1 = sdpvar(1);


S = sdpvar(2);
% Potential we found
% S = [1,0;0,0];

cons_SDP = [x1-xs g1/L]*S*[x1-xs g1/L].' - [x0-xs g0/L]*S*[x0-xs g0/L].';



cons_LIN = t1*(f1-fs) - t0*(f0-fs);


nu = sdpvar(nbPts,nbPts,'full');


for i = 1:nbPts
    for j = 1:nbPts
        if j ~= i  %&& (i~=2 || j~=3) 
            xi = XPEP(:,i); xj = XPEP(:,j);
            gi = GPEP(:,i); gj = GPEP(:,j);
            fi = FPEP(:,i); fj = FPEP(:,j);

            
            cons_SDP = cons_SDP - nu(i,j) * (gi*(xj-xi).' + 1/2/L*(gi-gj)*(gi-gj).');
            cons_LIN = cons_LIN - nu(i,j) * (fi - fj);
        end
    end
end


cons_SDP = (cons_SDP.' + cons_SDP)/2;
cons = cons_SDP <=0;
cons = cons+ (cons_LIN == 0);
cons = cons + (S >= 0);
cons = cons + (trace(S) <= 1);
cons = cons + (nu >= 0);



% maximize t1
obj =-t1;


solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
solverDetails=optimize(cons,obj,solver_opt);

"t1 = "+(double(t1))+", t0 + 2 lambda = "+(t0+2*lambda)