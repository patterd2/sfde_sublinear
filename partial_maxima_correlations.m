% Code to produce figure 6 in the paper of general sublinear SFDE results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function definitions
sigma = @(t,L_f,beta)  (L_f^(1/(1-beta)))...
    *((1-beta)*t)^(0.5*((1+beta)/(1-beta)))/sqrt(log(log(t+exp(1))));
f = @(t,beta) sign(t)*(abs(t))^beta; 
pow = @(x,alpha) sign(x).*abs(x).^alpha; % for powers of negative numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng = ('simdTwister');
% Input parameters for the FDE
beta = 0.5;
L_f = 1.05;
h = 0.01; % step size
X_0 = 1; % initial condition
epsilon = 0.07;

% set the terminal time of the simulation in "real time"
T = 100000;
% length of simulation in "discretised time" including initial interval
T_h = floor(T/h);
% create vector to store the solution
M_h = zeros(T_h+1,1);
X_h = zeros(T_h+1,1);
X_h(1,1) = X_0;
I_h = zeros(T_h+1,1);
W_n = sqrt(h)*randn(T_h,1); % normal increments for Brownian motion
% calculate the solution on [h, T_h] or (0,T]
for i = 1:T_h;
    I_h(i+1,1) = I_h(i,1) - h*(I_h(i,1) - f(X_h(i,1),beta));
    X_h(i+1,1) = X_h(i,1)+h*I_h(i,1)+sigma(i*h,L_f,beta)*W_n(i,1);
    M_h(i+1,1) = sigma(i*h,L_f,beta)*W_n(i,1);
end
% Plot s(nh)*X_h(nh)/Sigma(nh), with nonlinear transformation of the axes
t = 0:h:(T_h)*h;
s = (1+t).^(-3);
Sigma = (L_f*(1-beta)*t).^(1/(1-beta));
hold on;
plot(pow(t,epsilon),pow(transpose(s).*M_h./transpose(Sigma),epsilon),...
    'Color','g','LineWidth', 1.5);
plot(pow(t,epsilon),pow(transpose(s).*X_h./transpose(Sigma),epsilon),...
    'Color','r','LineWidth', 1.5);
hold on;
plot(pow(t,epsilon),pow(s*(L_f/(L_f-1)),epsilon),'Color','b','LineWidth',1.5);
hold on;
plot(pow(t,epsilon),pow(s*((L_f-1)/L_f),epsilon),'Color','m','LineWidth',1.5);
hold on;
plot(pow(t,epsilon),pow(-s*((L_f-1)/L_f),epsilon),'Color','m','LineWidth',1.5);
hold on;
plot(pow(t,epsilon),pow(-s*(L_f/(L_f-1)),epsilon),'Color','b','LineWidth',1.5);
hold on;
set(gca,'FontSize',22)
xlabel('Time','Interpreter','Latex')
set(gca,'YLim',[-1 1])
set(gca,'XLim',[1 max(pow(t,epsilon))]);
h = legend('$\frac{s(t)\,M(t)}{\Sigma(t)}$','$\frac{s(t)\,X(t)}{\Sigma(t)}$',...
    '$\pm \frac{L_f(\Sigma)}{L_f(\Sigma)-1}s(t)$','$\pm \frac{L_f(\Sigma)-1}{L_f(\Sigma)}s(t)$');
set(h,'Interpreter','latex')
positions = pow([0 10 50 250 1000 2500 10000 25000 100000],epsilon);
set(gca,'xTick',positions);
set(gca, 'xTickLabel',[0 10 50 250 1000 2500 10000 25000 100000]);
c = zeros(length(t),1);
hold on;
plot(t,c,'LineWidth',1,'Color','k');


