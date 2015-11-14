% Code to produce figure 6 in the paper of general sublinear SFDE results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function definitions
sigma = @(t,Lambda,beta)  (Lambda^(1/(1-beta)))...
    *((1-beta)*t)^(0.5*((1+beta)/(1-beta)))/sqrt(log(log(t+exp(1))));
f = @(t,beta) sign(t)*(abs(t))^beta; 
pow = @(x,alpha) sign(x).*abs(x).^alpha; % for powers of negative numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng = ('simdTwister');
% Input parameters for the FDE
beta = 0.5;
Lambda = 2;
h = 0.01; % step size
X_0 = 1; % initial condition

% set the terminal time of the simulation in "real time"
T = 1000;
% length of simulation in "discretised time" including initial interval
T_h = floor(T/h);
% create vector to store the solution

X_h = zeros(T_h+1,1);
X_h(1,1) = X_0;
I_h = zeros(T_h+1,1);
W_n = sqrt(h)*randn(T_h,1); % normal increments for Brownian motion
% calculate the solution on [h, T_h] or (0,T]
for i = 1:T_h;
    I_h(i+1,1) = I_h(i,1) - h*(I_h(i,1) - f(X_h(i,1),beta));
    X_h(i+1,1) = X_h(i,1)+h*I_h(i,1)+sigma(i*h,Lambda,beta)*W_n(i,1);
end
% Plot s(nh)*X_h(nh)/Sigma(nh), with nonlinear transformation of the axes
t = 0:h:(T_h)*h;
s = (1+t).^(-3);
Sigma = (Lambda*(1-beta)*t).^(1/(1-beta));
plot(pow(t,0.06),pow(transpose(s).*X_h./transpose(Sigma),0.06),...
    'Color','r','LineWidth', 1.5);
hold on;
plot(pow(t,0.06),pow(s*(Lambda/(Lambda-1)),0.06),'Color','b','LineWidth',1.5);
hold on;
plot(pow(t,0.06),pow(s*((Lambda-1)/Lambda),0.06),'Color','m','LineWidth',1.5);
hold on;
plot(pow(t,0.06),pow(-s*((Lambda-1)/Lambda),0.06),'Color','m','LineWidth',1.5);
hold on;
plot(pow(t,0.06),pow(-s*(Lambda/(Lambda-1)),0.06),'Color','b','LineWidth',1.5);
hold on;
set(gca,'FontSize',22)
xlabel('Time','Interpreter','Latex')
set(gca,'YLim',[-1 1])
set(gca,'XLim',[1 max(pow(t,0.06))]);
h = legend('$\frac{s(t)\,X_t(\omega)}{\Sigma(t)}$',...
    '$\pm \frac{\Lambda}{\Lambda-1}s(t)$','$\pm \frac{\Lambda-1}{\Lambda}s(t)$');
set(h,'Interpreter','latex')
positions = pow([0 10 50 250 1000 2500 10000 25000 100000],0.06);
set(gca,'xTick',positions);
set(gca, 'xTickLabel',[0 10 50 250 1000 2500 10000 25000 100000]);
c = zeros(length(t),1);
hold on;
plot(t,c,'LineWidth',1,'Color','k');


