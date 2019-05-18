function result = ADMM(n, R, rho, Z_0, Phi_0, alpha, beta, epsilon)
%% Alternating Direction Method of Multipliers (ADMM) 
% for a Special Class of Bilinear Programming:
%    min     2*tr(X^TR) + tr(Z^TXR)
%    s.t.    Xe = rho, tr(X)=0, 
%            Z^Te = rho, Z>=0,
%            X = Z.

% Input:
% n: Dimension of Problem
% R: Symmetric Constant Matrix (Distance??), with diagonal entries 0
%    and all entries non-negative
% rho: Constant Vector, with all components positive
% Z_0: Initial Splitting Matrix
% Phi_0: Initial Multiplier for Constraint 'X = Z'
% alpha: Over-relaxation Paramater (0~1)
% beta: Penalty Parameter, >0 (say, 1e3)
% epsilon: Tolerance

% Output: a struct 'result'
% _.k: Iteration Number
% _.X: Output Transport Map X
% _.Z: Output Splitting Variable Z
% _.f: Objective Value
% _.t: Running Time (s)
% _.E: KKT Violation

%% Initialization
e = ones(n, 1);
Phi = Phi_0;
Z = Z_0;
T = 1; S = 1;% T -- primal residual, S -- dual residual
E = T + S; % KKT Violation
k = 0; % iteration number
options = optimoptions('quadprog','Display','off',...
    'ConstraintTolerance',1e-10);

%% ADMM Scheme
tic
while E > epsilon
    %% X subproblem
    % Problem Description:
    %   min  2tr(X^TR) + tr(Z^TXR) - tr(Phi'(X-Z)) + beta/2||X - Z||_F^2
    %   s.t. Xe = rho, tr(X)=0
    
    temp_1 = R*e; temp_2 = beta*rho;
    M_1 = 2*temp_1 + Z*temp_1 - Phi*e - beta*Z*e + temp_2;
    m_1 = tr(Z, R) - trace(Phi) - beta*trace(Z);
    mu = 1/(n-1)*(-1/n*e'*M_1 + m_1);
    % multiplier for constraint 'tr(X)=0'
    lambda_1 = 1/n*(M_1 - e*mu);
    % multiplier for constraint 'Xe = rho'
    X = -1/beta*(2*R + Z*R - lambda_1*e' - mu*eye(n) - Phi - beta*Z);

    %% Z subproblem
    % Problem Description
    %   min  tr(Z^TXR) - tr(Phi'(X-Z)) + beta/2||X - Z||_F^2
    %   s.t. Z^Te = rho, Z>=0
    
    lambda_2 = 1/n*(R*X'*e + Phi'*e - beta*X'*e + temp_2);
    % multiplier for constraint 'Z^Te = rho'
    tildeZ = -1/beta*(X*R + Phi - beta*X - e*lambda_2');
    % solution to equality constrained problem, ignoring Z >= 0
    
    % Column Subproblem
    % Problem Description:
    % min   \sum_{i=1}^nt_i^2
    % s.t.  t_i>=-z_i,i=1,...,n
    %       \sum_{i=1}^nt_i=0
    tildeZ = Column_Alg(tildeZ, n);
    
    %% Update Multipliers
    Phi = Phi - beta*alpha*(X - tildeZ); % over-relaxation update for multiplier
    
    %% Draw
%     if k==0
%         f = 2*tr(X,R)+trace(tildeZ'*X*R);
%         T = norm(X-tildeZ,inf);
%         S = norm(X-tildeZ,inf);
%         figure(1)
%         plot(0,log10(T),'ro-'),hold on
%         figure(2)
%         plot(0,log10(S),'bo-'),hold on
%         figure(3)
%         plot(0,f,'go-'),hold on
%     elseif k > 0
%         T_old = T; S_old = S; f_old = f;
%         f = 2*tr(X,R)+trace(tildeZ'*X*R);
%         T = norm(X-tildeZ,inf);
%         S = norm((tildeZ-Z)*(beta*eye(n)-R),inf);
%         figure(1)
%         plot([k-1,k],[log10(T_old),log10(T)],'ro-'),hold on
%         figure(2)
%         plot([k-1,k],[log10(S_old),log10(S)],'bo-'),hold on
%         figure(3)
%         plot([k-1,k],[f_old,f],'go-'),hold on
%     end
    %% Update residual
    %f = 2*trace(X'*R) + trace(tildeZ'*X*R);
    T = norm(X - tildeZ, inf);
    S = norm((tildeZ - Z)*(beta*eye(n)-R), inf);
    E = T + S;
    k = k + 1;
    Z = tildeZ;
end
t = toc;

% Pack up

f = 2*trace(X'*R) + trace(Z'*X*R);
result = struct('k',k,'X',X,'Z',Z,'f',f,'t',t,'E',E);

end

function f = tr(A, B)
% Calculate the trace of the product of two matrices

n = size(A, 1);
f = 0;
for i = 1 : n
    f = f + A(i, :)*B(:, i);
end
end

function tildeX = Column_Alg(X, n)
% Column algorithm for subproblem

tildeX = X;
for j = 1 : n
    if min(tildeX(:, j)) >= 0
        continue;
    else
        ind = tildeX(:, j) <= 0;
        ind_m = tildeX(:, j) > 0;
        t = zeros(n, 1);
        t = max(t, -tildeX(:, j));
        while sum(ind) < n
            bar_t = -sum(t)/sum(ind_m);
            if tildeX(ind_m, j) >= -bar_t
                t(ind_m) = bar_t;
                tildeX(:, j) = tildeX(:, j) + t;
                break;
            else
                ind_a = tildeX(:, j) < -bar_t;
                ind_a = ind_a - ind;
                ind_a = logical(max(0, ind_a));
                t(ind_a) = -tildeX(ind_a, j);
                ind = logical(ind + ind_a);
                ind_m = logical(ind_m - ind_a);
            end
        end
    end
end
end