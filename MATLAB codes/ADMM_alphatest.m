for n=[3,4,5]
    epsilon = 1e-5;
    beta = 1e2;
    
    Rname = strcat('R_', num2str(n));
    rhoname = strcat('rho_', num2str(n));
    Phiname = strcat('Phi_', num2str(n));
    Zname = strcat('Z_', num2str(n));
    load(Rname)
    load(rhoname)
    load(Phiname)
    load(Zname)
    
    for alpha = 0.1 : .1 : 1.0
        result = ADMM(n, R, rho, Z_0, Phi_0, alpha, beta, epsilon);
        save(strcat('alpha', num2str(10*alpha), 'n', num2str(n)), 'result')
    end
    fprintf('n = %d complete\n', n)
end

for n=[20,30,40]
    epsilon = 1e-4;
    beta = 1e4;
    
    Rname = strcat('R_', num2str(n));
    rhoname = strcat('rho_', num2str(n));
    Phiname = strcat('Phi_', num2str(n));
    Zname = strcat('Z_', num2str(n));
    load(Rname)
    load(rhoname)
    load(Phiname)
    load(Zname)
    
    for alpha = 0.1 : .1 : 1.0
        result = ADMM(n, R, rho, Z_0, Phi_0, alpha, beta, epsilon);
        save(strcat('alpha', num2str(10*alpha), 'n', num2str(n)), 'result')
    end
    fprintf('n = %d complete\n', n)
end