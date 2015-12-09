% Code to produce figure 5 in the paper of general sublinear SFDE results
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
L_f = [0.5 1.5 5];
h = 0.01; % step size
X_0 = 1; % initial condition, not important, choose nonzero
line_colours = char('r', 'b', 'm');

% set the terminal time of the simulation in "real time"
T = 100000;
% length of simulation in "discretised time" including initial interval
T_h = floor(T/h);
% create vector to store the solution

for j = 1:length(L_f);
    X_h = zeros(T_h+1,1);
    X_h(1,1) = X_0;
    I_h = zeros(T_h+1,1);
    W_n = sqrt(h)*randn(T_h,1); % normal increments for Brownian motion
    
    % calculate the solution on [h, T_h] or (0,T]
    for i = 1:T_h;
        I_h(i+1,1) = I_h(i,1) - h*(I_h(i,1) - f(X_h(i,1),beta));
        X_h(i+1,1) = X_h(i,1)+h*I_h(i,1)+...
            sigma(i*h,L_f(1,j),beta)*W_n(i,1);
    end
    
    t = 0:h:(T_h)*h;
    % Plot F(|X_h(nh)|)/nh
    plot(t,(((abs(X_h)).^(1-beta))/(1-beta))./transpose(t),'Color',...
        line_colours(j,1),'LineWidth',1.5);
    hold on;
end
hold on;
set(gca,'FontSize',22)
xlabel('Time','Interpreter','Latex')
set(gca,'XLim',[3 T])
set(gca,'YLim',[0 8])
set(gca,'XTick',0:T/4:T)
h=legend('$L_f(\Sigma)=0.5$','$L_f(\Sigma)=1.5$','$L_f(\Sigma)=5$');
set(h,'Interpreter','latex');



