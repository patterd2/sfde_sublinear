% Code to produce figure 3 in the paper of general sublinear SFDE results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define functions
sigma = @(t,alpha) sqrt(alpha+1)*t^(alpha/2); % diffusion coefficient
f = @(t,beta) sign(t)*(abs(t))^beta; % nonlinearity
pow = @(x,alpha) sign(x).*abs(x).^alpha; % for powers of negative numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng = ('simdTwister');
% Input parameters for the FDE
beta = 0.5;
h = 0.01; % step size
X_0 = 1; % initial condition
alpha = 4;
epsilon = 0.07;

% set the terminal time of the simulation in "real time"
T = 250000;
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
    X_h(i+1,1) = X_h(i,1)+h*I_h(i,1)+sigma(i*h,alpha)*W_n(i,1);
end
% Plot s(nh)*X_h(nh)/Sigma(nh), with nonlinear transformation of the axes
t = 0:h:(T_h)*h;
s = (1+t).^(-3);
Sigma = sqrt(2*(t.^(alpha+1)).*log(log(t.^(alpha+1)+exp(1))));
plot(pow(t,epsilon),pow(transpose(s).*X_h./transpose(Sigma),epsilon),...
    'Color','r','LineWidth', 1.5);
hold on;
plot(pow(t,epsilon),pow(s,epsilon),'Color','b','LineWidth',1.5);
plot(pow(t,epsilon),pow(-s,epsilon),'Color','b','LineWidth',1.5);
set(gca,'FontSize',22)
xlabel('Time','Interpreter','Latex')
set(gca,'YLim',[-1 1])
set(gca,'XLim',[1 max(pow(t,epsilon))]);
h = legend('$\frac{s(t)\,X(t)}{\Sigma(t)}$','$\pm s(t)$');
set(h,'Interpreter','latex');
positions = pow([0 10 50 250 1000 3000 10000 30000 100000 250000],epsilon);
set(gca,'xTick',positions);
set(gca, 'xTickLabel',[0 10 50 250 1000 3000 10000 30000 100000 250000]);
c = zeros(length(t),1);
plot(t,c,'LineWidth',1,'Color','k');


