options = optimoptions('fmincon','Algorithm','sqp','ConstraintTolerance',1e-8);
for n=[3,4,5]
    Rname = strcat('R_', num2str(n));
    rhoname = strcat('rho_', num2str(n));
    Phiname = strcat('Phi_', num2str(n));
    Zname = strcat('Z_', num2str(n));
    load(Rname)
    load(rhoname)
    load(Phiname)
    load(Zname)
    I = eye(n);
    e = ones(n,1);
    tic
    [xy,fval] = fmincon(@(x) myfun(x,R),Z_0(:),[],[],...
        [kron(e',I);kron(I,e');I(:)'],[rho;rho;0],zeros(n^2,1),...
        [],[],options);
    toc
    fprintf('%.4f\n',fval)
end
    
options = optimoptions('fmincon','Algorithm','sqp','ConstraintTolerance',1e-6);
for n=[20,30,40,50,60,70,80]
    Rname = strcat('R_', num2str(n));
    rhoname = strcat('rho_', num2str(n));
    Phiname = strcat('Phi_', num2str(n));
    Zname = strcat('Z_', num2str(n));
    load(Rname)
    load(rhoname)
    load(Phiname)
    load(Zname)
    I = eye(n);
    e = ones(n,1);
    tic
    [xy,fval] = fmincon(@(x) myfun(x,R),Z_0(:),[],[],...
        [kron(e',I);kron(I,e');I(:)'],[rho;rho;0],zeros(n^2,1),...
        [],[],options);
    toc
    fprintf('%.4f\n',fval)
end

function f = myfun(x,R)

n = size(R,1);
X = reshape(x,n,n);
f = 2*trace(X*R) + trace(X'*X*R);
end