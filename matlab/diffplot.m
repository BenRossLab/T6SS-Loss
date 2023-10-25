function diffplot(tt,pp,rr,ss,fr)

% Takes genotype abundances through time from spatial simultaion
% Generates plot of spatial distribution of genotype abundances at frame fr

% Colors associated with each phenotype
CP = [0.8470588235294118, 0.10588235294117647, 0.3764705882352941];     % Producer: Magenta
CR = [0.11764705882352941, 0.5333333333333333, 0.8980392156862745];     % Resistant: Blue
CS = [0.8823529411764706, 0.7568627450980392, 0.027450980392156862];    % Sensitive: Yellow

% Space intervals
N=51;   % length of grid
dx = 1/(N-1);
x = [0:dx:1]; y = [0:dx:1];

% Calculates RGB matrix to plot
pp=reshape(pp(fr+1,:),N,N);
ss=reshape(ss(fr+1,:),N,N);
rr=reshape(rr(fr+1,:),N,N);
C(:,:,1)=min(1,pp*CP(1)+rr*CR(1)+ss*CS(1));
C(:,:,2)=min(1,pp*CP(2)+rr*CR(2)+ss*CS(2));
C(:,:,3)=min(1,pp*CP(3)+rr*CR(3)+ss*CS(3));

% Plots abundances
image(x,y,C);
axis off
text(0.4,0.1,['Time ' num2str(tt(fr+1))],'Color','w','FontSize',15,'FontWeight','bold')


