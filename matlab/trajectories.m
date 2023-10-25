function xx=trajectories(ap,ar,g,x0)

% Takes parameters for cost of production (ap), cost of resistance (ar), and density (g) 
% Starts at initial abundances (x0)
% Integrates differential equations and plots results
% Returns species abundances through time

% Colors associated with each phenotype
CP = [0.8470588235294118, 0.10588235294117647, 0.3764705882352941];     % Producer: Magenta
CR = [0.11764705882352941, 0.5333333333333333, 0.8980392156862745];     % Resistant: Blue
CS = [0.8823529411764706, 0.7568627450980392, 0.027450980392156862];    % Sensitive: Yellow

% Integrates differential equations
p=[ap ar g];
opts = odeset('Events',@(t,x) events(t,x,p));
[tt,xx]=ode113(@(t,x) equations(t,x,p),[0 300],x0,opts);
tt=tt/log(2);

% Plots results
hold on
plot(tt,xx(:,1),'Color',CP,'LineWidth',3)
plot(tt,xx(:,2),'Color',CR,'LineWidth',3)
plot(tt,xx(:,3),'Color',CS,'LineWidth',3)
xlabel('Time (doublings)','FontSize',20)
ylabel('Fraction of population','FontSize',20)
legend('Producer','Resistant','Sensitive','FontSize',15)
axis([0 tt(end) 0 1.1])
set(gca,'FontSize',15)

% Differential equations
function dx=equations(t,x,p)
    X=@(x) x(1)*(1-p(1))+x(2)*(1-p(2))+x(3);
    dx=zeros(3,1);
    dx(1)=x(1).*(1-p(1)-X(x));
    dx(2)=x(2).*(1-p(2)-X(x));
    dx(3)=x(3).*(1-X(x)-p(3)*x(1));
end

% Stops simulation if steady-state is reached
function [value,isterminal,direction]=events(t,x,p)
    dx=equations(t,x,p);
    value=norm(dx./x)-1e-3;
    isterminal=1;
    direction=0;
end
end