% Code to produce figure 2 in the paper of general sublinear SFDE results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define functions
sigma = @(t) exp(-t); % diffusion coefficient, in L^2
f = @(t,beta) sign(t)*(abs(t))^beta; % nonlinearity
pow = @(x,alpha) sign(x).*abs(x).^alpha; % for powers of negative numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng = ('simdTwister');
% Input parameters for the FDE
beta = [0.25 0.5 0.75];
h = 0.01; % step size
X_0 = 1; % initial condition
line_colours = char('r', 'b', 'm');

% set the terminal time of the simulation in "real time"
T = 10000;
% length of simulation in "discretised time" including initial interval
T_h = floor(T/h);
% create vector to store the solution

for j = 1:3;
    X_h = zeros(T_h+1,1);
    X_h(1,1) = X_0;
    I_h = zeros(T_h+1,1);
    W_n = sqrt(h)*randn(T_h,1); % normal increments for Brownian motion
    
    % calculate the solution on [h, T_h] or (0,T]
    for i = 1:T_h;
        I_h(i+1,1) = I_h(i,1) - h*(I_h(i,1) - f(X_h(i,1),beta(1,j)));
        X_h(i+1,1) = X_h(i,1)+h*I_h(i,1)+sigma(i*h)*W_n(i,1);
    end
    
    t = 0:h:(T_h)*h;
    % Plot F(|X_h(nh)|)/nh
    plot(t,(((abs(X_h)).^(1-beta(1,j)))/(1-beta(1,j)))./transpose(t),...
        'Color',line_colours(j,1),'LineWidth',1.5);
    hold on;
end
set(gca,'FontSize',22)
xlabel('Time','Interpreter','Latex')
set(gca,'XLim',[3 T])
set(gca,'YLim',[0.8 1.2])
set(gca,'YTick',[0.8 1 1.2])
set(gca,'XTick',0:T/4:T)
c = ones(length(t),1);
plot(t,c,'LineWidth',1,'Color','k');
h = legend('$\beta=0.25$','$\beta=0.5$','$\beta=0.75$','Location','northeast');
set(h,'Interpreter','latex');



