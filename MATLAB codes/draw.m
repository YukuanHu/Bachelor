for n=[3,4,5,20,30,40]
    for i = 1:10
        load(strcat('alpha',num2str(i),'n',num2str(n)))
        fig = imagesc(result.X);
        alphaname = strcat('$\alpha=',num2str(i/10),'$');
        title(alphaname,'Interpreter','latex')
        saveas(fig,strcat('C:\Users\hyk\Desktop\同济大学毕业论文模板\figures\alpha',num2str(i),'n',num2str(n),'.jpg'))
        close all
    end
end