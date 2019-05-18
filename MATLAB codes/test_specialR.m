p = 3; q = 4; % p \ne q; p, q < n
alpha = 1;
beta = 1e3;
epsilon = 1e-8;

for n = 5 : 5 : 20
    R = zeros(n, n);
    R(p, q) = 1; R(q, p) = 1;
    rho = ones(n, 1);
    Z_0 = rand(n, n); 
    Phi_0 = rand(n, n); 
    save(strcat('R_',num2str(n),'special'),'R')
    save(strcat('rho_',num2str(n),'special'),'rho')
    save(strcat('Z_',num2str(n),'special'),'Z_0')
    save(strcat('Phi_',num2str(n),'special'),'Phi_0')
    result = ADMM(n, R, rho, Z_0, Phi_0, alpha, beta, epsilon);
    fprintf('迭代数: %d\n所耗时间 (s): %.4f\nKKT违反度: %.2e\n目标值: %.2e\n',...
        result.k,result.t,result.E,result.f)
end