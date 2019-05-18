for n = [3,4,5,20,30,40,50,60,70,80]
    R = abs(randn(n));
    for i = 1:n
        R(i,i)=0;
    end
    R = (R+R')/2;
    rho = abs(randn(n,1));
    Phi_0 = randn(n);
    Z_0 = randn(n);
    save(strcat('R_',num2str(n)),'R')
    save(strcat('rho_',num2str(n)),'rho')
    save(strcat('Phi_',num2str(n)),'Phi_0')
    save(strcat('Z_',num2str(n)),'Z_0')
end
