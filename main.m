
%% choose the correct flags
linear_variation=1; % equal theta
Random_strengths=0; % unequal theta

% More flags
deco=1; % 0- no decoherence, 1-decoherence
other_rates=0; % good to keep 1 for high temp
depolarization=1; % 0 depolarization channel is off, 1 on

%% Other User parameters

    
nbs=26; % No. of beam splitters
last=180;  % No. of experiments (M)

%% Simulation
Arbitrary_N_theta_sim

%% Experimental data
load('net_p0.mat')
load('eta_exp.mat')
load('pr_exp.mat')
load('nr_exp.mat')

figure
tiledlayout (2,2)

nexttile %1
   Y1=linspace(0,1,n_experiments);
   X1=linspace(1,n_bombs,n_bombs);
   [X,Y]=meshgrid(X1,Y1);

   
   surf(X,Y,eta_exp)
   view ([90, -90])
   shading 'flat'
   colorbar
   caxis ([0 1])
   ylim ([0 1])
   xlim ([1 n_bombs])
    title("\eta_c Experiment")
    xlabel('N');
    ylabel('\theta/\pi');
%%%%%
   nexttile %2
   Y1=linspace(0,1,n_experiments);
   X1=linspace(1,n_bombs,n_bombs);
   [X,Y]=meshgrid(X1,Y1);
   
   surf(X,Y,abs(net_p0))
   view ([90, -90])
   shading 'flat'
   colorbar
   caxis ([0 1])
   ylim ([0 1])
   xlim ([1 n_bombs])

    title("p_0 Experiment")
    xlabel('N');
    ylabel('\theta/\pi');
  %%%%% %3
   nexttile
   Y1=linspace(0,1,n_experiments);
   X1=linspace(1,n_bombs,n_bombs);
   [X,Y]=meshgrid(X1,Y1);
   
   surf(X,Y,pr_exp)
   view ([90, -90])
   shading 'flat'
   colorbar
   caxis ([0 1])
   ylim ([0 1])
   xlim ([1 n_bombs])

    title("PR(\theta) Experiment")
    xlabel('N');
    ylabel('\theta/\pi');

    %%%%% %3
   nexttile
   Y1=linspace(0,1,n_experiments);
   X1=linspace(1,n_bombs,n_bombs);
   [X,Y]=meshgrid(X1,Y1);
   
   surf(X,Y,nr_exp)
   view ([90, -90])
   shading 'flat'
   colorbar
   caxis ([0 1])
   ylim ([0 1])
   xlim ([1 n_bombs])

    title("NR(\theta) Experiment")
    xlabel('N');
    ylabel('\theta/\pi');

    %%
    figure 
    tiledlayout(2,1)

    nexttile
    plot(X1, mean(net_p0), X1, net_p0(1,:), X1, net_p0(end,:))
 title("Experiment: E[p_0] (blue), p_0(\theta=0) (red) and p_0(\theta=\pi) (yellow)" )
    xlabel('N');
    ylabel('D_0 outcome');

 nexttile
  plot(X1,std(net_p0))
title("Standard deviation" )
    xlabel('N');
    ylabel('\sigma[p_0]');