% Code to analse the characteristic equation for RV results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 1.5;
betas = 0.5:0.0005:0.9999;
N = length(betas);
etas = ones(1,N);
for i = 1:N;
    beta = betas(1,i);
    f = @(x) x - (x.^beta)/lambda - 1;
    etas(1,i) = fsolve(f,1);
end
plot(betas,etas,'LineWidth',2);
set(gca,'FontSize',22)
hold on;
c = (lambda/(lambda-1))*ones(1,N);
plot(betas,c,'LineWidth',2);
xlabel('$\beta$','Interpreter','Latex')
set(gca,'xTick',[0.5 0.6 0.7 0.8 0.9 1]);
set(gca, 'xTickLabel',[0.5 0.6 0.7 0.8 0.9 1]);
h = legend({'$\eta(\beta)$','$\frac{\lambda}{\lambda-1}$'},'Location','northwest');
set(h,'Interpreter','latex')


