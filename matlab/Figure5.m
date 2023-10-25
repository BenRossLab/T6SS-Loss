% Plots figure 5 as in paper
% Uses corresponding spatial simulation results traj1, traj2, traj2b, and traj3

% Parameters
ap=0.05;            % Cost of production
ar=0.02;            % Cost of resistance
gg=[0.1 0.2 0.5];   % Density = killing*quality/growth

x0=[0.3 0.05 0.65]; % Initial conditions (producer, resistant, sensitive)

much={'Low','Moderate','High'};

figure
for i=1:3

    g=gg(i);
    
    % Plots deterministic trajectories
    subplot(2,3,3+i)
    xx=trajectories(ap,ar,g,x0);

    % Plots corresponding triangle plots
    subplot(2,3,i)
    triplot(ap,ar,g,xx)

    title([much{i} ' Density'],'FontSize',30)
    
end

figure
names={'traj1','traj2b','traj2','traj3'};
for i=1:4
    
    load(names{i})
    
    % Plots spatial simulations
    subplot(4,4,i)
    diffplot(tt,pp,rr,ss,10)
    subplot(4,4,4+i)
    diffplot(tt,pp,rr,ss,100)
    subplot(4,4,8+i)
    diffplot(tt,pp,rr,ss,300)
    subplot(4,4,12+i)
    diffplot(tt,pp,rr,ss,600)

end