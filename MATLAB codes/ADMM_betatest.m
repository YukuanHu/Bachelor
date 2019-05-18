for n=[3,4,5]
    epsilon = 1e-8;
    alpha = 1;
    
    Rname = strcat('R_', num2str(n));
    rhoname = strcat('rho_', num2str(n));
    Phiname = strcat('Phi_', num2str(n));
    Zname = strcat('Z_', num2str(n));
    load(Rname)
    load(rhoname)
    load(Phiname)
    load(Zname)
    
    for beta = [1e2,1e3,1e4,1e5]
        result = ADMM(n, R, rho, Z_0, Phi_0, alpha, beta, epsilon);
        save(strcat('beta', num2str(beta), 'n', num2str(n)), 'result')
    end
    fprintf('n = %d complete\n', n)
end

for n=[20]
    epsilon = 1e-6;
    alpha = 1;
    
    Rname = strcat('R_', num2str(n));
    rhoname = strcat('rho_', num2str(n));
    Phiname = strcat('Phi_', num2str(n));
    Zname = strcat('Z_', num2str(n));
    load(Rname)
    load(rhoname)
    load(Phiname)
    load(Zname)

    for beta = [1e3,1e4,1e5]
        result = ADMM(n, R, rho, Z_0, Phi_0, alpha, beta, ...
            epsilon);
        save(strcat('beta', num2str(beta), 'n', num2str(n)), 'result')
    end
    fprintf('n = %d complete\n', n)
end

for n=[30]
    epsilon = 1e-6;
    alpha = 1;
    
    Rname = strcat('R_', num2str(n));
    rhoname = strcat('rho_', num2str(n));
    Phiname = strcat('Phi_', num2str(n));
    Zname = strcat('Z_', num2str(n));
    load(Rname)
    load(rhoname)
    load(Phiname)
    load(Zname)
    % n=30 beta=1e3ª·—≠ª∑
    for beta = [1e4,1e5]
        result = ADMM(n, R, rho, Z_0, Phi_0, alpha, beta, ...
            epsilon);
        save(strcat('beta', num2str(beta), 'n', num2str(n)), 'result')
    end
    fprintf('n = %d complete\n', n)
end