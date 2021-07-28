% Code to reproduce Figure 3.2
clc; clear all;

verbose   = 1;
tolerance = 1e-8;


L=1;

lambda = 1.5/L;
Ns     = 1:10;
res    = zeros(length(Ns),1);
k = 1;
for N = Ns


    dimG = N+3;dimF = N+2;

    xs = zeros(dimG,1); xs(1,1) = 1;
    x  = zeros(dimG,N+1);x(2,1) = 0;
    gs = zeros(dimG,1);
    g  = zeros(dimG,N+1);g(3:N+3,:) = eye(N+1);

    %gradient updates
    for i = 1:N
       x(:,i+1) = x(:,i) - lambda*g(:,i); 
    end

    fs = zeros(dimF,1);fs(1,1) = 1;
    f  = zeros(dimF,N+1);f(2:N+2,:) = eye(N+1);

    G = sdpvar(dimG);
    F = sdpvar(1,dimF);

    obj = F*(f(:,N+1)-fs);

    %initial condition
    cons = (x(:,1)-xs).'*G*(x(:,1)-xs) <= 1;

    nbPts = N+2;
    XPEP  = [xs x];
    GPEP  = [gs g];
    FPEP  = [fs f];

    for i = 1:nbPts 
        for j = 1:nbPts
            if i~=j
                xi = XPEP(:,i); xj = XPEP(:,j);
                gi = GPEP(:,i); gj = GPEP(:,j);
                fi = FPEP(:,i); fj = FPEP(:,j);
                cons = cons + (F*(fi-fj) + gi.'*G*(xj-xi) + 1/2/L*(gi-gj).'*G*(gi-gj) <=0);
            end
        end
    end

    cons = cons + (G >= 0);

    solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
    solverDetails=optimize(cons,-obj,solver_opt);
    
    res(k) = double(obj);
    k      = k+1;
end

plot(Ns,res);
hold on
plot(Ns,L/2*max((1-L*lambda).^(2*Ns),1./(2*Ns*L*lambda+1)),'LineStyle','--');
legend(["C(N,1)", "conjuecture Drori and Teboulle"]);