%clear all
%clf
clear i
clear sim_plot_data pop0 pop1 pop2 rhof
%%
GAMMA10 = 1*1/3.4*1e6; % Relaxation rate for 0-1 transition
GAMMA10phi = 0.178e6;  % Pure dephasing rate for 0-1 transition
GAMMA21 = 1.152e6;     % Relaxation rate for 1-2 transition
GAMMA21phi = 1.8259e6; % Pure dephasing rate for 1-2 transition
GAMMA20phi = 1.7e6;    % Pure dephasing rate for 0-2 transition
gamma10=1/2*GAMMA10 + GAMMA10phi;   % Dephasing rate for the 0-1 transition
gamma20=1/2*GAMMA21 + GAMMA20phi;   % Dephasing rate for the 0-2 transition
gamma21=1/2*(GAMMA21 + GAMMA10) + GAMMA21phi;   % Dephasing rate for the 1-2 transition
%%
GAMMA10/2/pi*1e-3
GAMMA21/2/pi*1e-3
GAMMA10phi/2/pi*1e-3
GAMMA21phi/2/pi*1e-3
GAMMA20phi/2/pi*1e-3
%%
    psi_i = [1;0;0];
    rho0 = psi_i*psi_i';
    nt=200;
    fraction=1/nt;
    duration=56e-9;
    dt=fraction*duration;
    
    
%     nbs=26; %51;
%     last=180;  % linear
n_bombs = nbs-1;
n_experiments=last;

% % choose the correct flags
% linear_variation=1; % equal theta
% Random_strengths=0; % unequal theta

if (linear_variation==1)
strength_range = linspace(0,1.0*pi,n_experiments+1);
elseif (Random_strengths==1)
strength_range = rand(n_experiments+1,n_bombs)*pi*1.0;
end
%% More flags
%deco=1;
%other_rates=0; % good to keep 1 for high temp
temp_change=0; % always keep at 0
%depolarization=1; % 0 depolarization cahnnel is off, 1 on
ideal_case=0;

  hb = 6.62607004*10^(-34);
  K = 1.380649*10^(-23);
  v01=7.2e9;%5*10^9;
  v12=6.85e9;%4.65*10^9;
  v02=v01+v12;
  initial_temperature=0.075;
  T=initial_temperature*(1-ideal_case)+0.001*ideal_case;%0.098;%0.098;
  pi1=exp(hb*v01/K/T)/(1+exp(hb*v01/K/T)+exp(-hb*v12/K/T));
  pi2=exp(-hb*v01/K/T)/(1+exp(-hb*v01/K/T)+exp(-hb*v02/K/T));
  pi3=exp(-hb*v02/K/T)/(1+exp(-hb*v01/K/T)+exp(-hb*v02/K/T));
  
  temperature_range = T+(temp_change)*linspace(1,100,n_experiments+1).*linspace(1,1,n_bombs)'*1e-3;
  rho0 = [pi1,0,0; 0,pi2,0; 0,0,pi3]
  
  dpol=0*depolarization*0.0012*fraction;
  dpol_mat=linspace(1/2,1/26,n_bombs)*2*dpol;
  dpol2_mat = linspace(0,2*2,length(strength_range))*0.0009*fraction*depolarization;
%%
avg_number=1;%10
sim_plot_data=zeros(length(strength_range),n_bombs,avg_number);
avg_counter=0;
for data_ind=1:avg_number
    
    for ind3=1:length(strength_range)
        
        %dpol2=dpol2_mat(ind3);

        for ind1=1:n_bombs
            final_time=0;
             angle = pi/(ind1+1)*fraction*1; %0.996; %0.995;            
            
   %%%% new addition on decoherence %%%%
  Temperature = temperature_range(ind1,ind3);
  relaxation_temp_01 = exp(-hb*v01/K/Temperature);
  relaxation_temp_12 = exp(-hb*v12/K/Temperature);
  relaxation_temp_02 = exp(-hb*v02/K/Temperature);
  %%%%%%%%%%%%% other rates, significant only at high temperatures %%%%%
GAMMA01 = relaxation_temp_01*GAMMA10;
GAMMA01phi = 0*GAMMA10phi;
GAMMA12 = relaxation_temp_12*GAMMA21;
GAMMA12phi = 0*GAMMA21phi;
GAMMA02phi = 0*GAMMA20phi;
gamma01 = 1/2*(GAMMA01 + GAMMA12) + GAMMA01phi;% + 1/4*(GAMMA12phi + GAMMA02phi);   % Dephasing rate for the 0-1 transition
gamma02 = 1/2*GAMMA01 + GAMMA02phi;% + 1/4*(GAMMA12phi + GAMMA01phi);  % Dephasing rate for the 0-2 transition
gamma12 = 1/2*GAMMA12 + GAMMA12phi;% + 1/4*(GAMMA02phi + GAMMA01phi);  % Dephasing rate for the 1-2 transition
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        rho=rho0;
        for ind_t=1:nt
            bs = [cos(angle/2), -sin(angle/2), 0; sin(angle/2), cos(angle/2), 0; 0,0,1];
             
            rho1 = bs'*rho*bs;
            Lrho = [GAMMA10*rho1(2,2), -gamma10*rho1(1,2), -gamma20*rho1(1,3);
                -gamma10*rho1(2,1), GAMMA21*rho1(3,3)-GAMMA10*rho1(2,2),-gamma21*rho1(2,3);
                -gamma20*rho1(3,1), -gamma21*rho1(3,2), -GAMMA21*rho1(3,3)]; %Lindblad term
            
            LrhoNEW = [-GAMMA01*rho1(1,1), -gamma01*rho1(1,2), -gamma02*rho1(1,3);
                -gamma01*rho1(2,1), -GAMMA12*rho1(2,2)+GAMMA01*rho1(1,1),-gamma12*rho1(2,3);
                -gamma02*rho1(3,1), -gamma12*rho1(3,2), GAMMA12*rho1(2,2)] ; %Lindblad term
            
            rho = (1-dpol_mat(ind1))*(rho1 + deco*(Lrho+other_rates*LrhoNEW)*dt) + dpol_mat(ind1)*eye(3)/3;   % addition of the rotated density matrix with Lindblad (dissipative) term
            final_time=final_time+dt;
        end

            for ind2=1:ind1
                if (Random_strengths==1)
                    strength_rand= strength_log(ind1,ind2,ind3,data_ind); % rand;
                     strength=strength_rand*pi*fraction*1.1;
                    strength_log(ind1,ind2,ind3,data_ind)=strength_rand;
                  %  
                 dpol2=strength_rand*4*0.9e-3*fraction*depolarization;
                elseif  (linear_variation==1)
                    strength = strength_range(ind3)*fraction*1.0;
                    dpol2=dpol2_mat(ind3);
                end
                
                for ind_t=1:nt
                    bomb = [1,0,0;0,cos(strength/2),-sin(strength/2); 0, sin(strength/2), cos(strength/2)];
                    rho1 = bomb'*rho*bomb;
                    Lrho = [GAMMA10*rho1(2,2), -gamma10*rho1(1,2), -gamma20*rho1(1,3);
                        -gamma10*rho1(2,1), GAMMA21*rho1(3,3)-GAMMA10*rho1(2,2),-gamma21*rho1(2,3);
                        -gamma20*rho1(3,1), -gamma21*rho1(3,2), -GAMMA21*rho1(3,3)]; %Lindblad term
                    
                    LrhoNEW = [-GAMMA01*rho1(1,1), -gamma01*rho1(1,2), -gamma02*rho1(1,3);
                        -gamma01*rho1(2,1), -GAMMA12*rho1(2,2)+GAMMA01*rho1(1,1),-gamma12*rho1(2,3);
                        -gamma02*rho1(3,1), -gamma12*rho1(3,2), GAMMA12*rho1(2,2)] ; %Lindblad term
                    
                    rho = (1-2*dpol2)*(rho1 + 2*deco*(Lrho+other_rates*LrhoNEW)*dt) + 2*dpol2*eye(3)/3;   % addition of the rotated density matrix with Lindblad (dissipative) term
                final_time=final_time+2*dt;
                end
                
                for ind_t=1:nt
                    angle2=angle*1;
                    bs = [cos(angle2/2), -sin(angle2/2), 0; sin(angle2/2), cos(angle2/2), 0; 0,0,1];
                  
                    rho1 = bs'*rho*bs;
                    Lrho = [GAMMA10*rho1(2,2), -gamma10*rho1(1,2), -gamma20*rho1(1,3);
                        -gamma10*rho1(2,1), GAMMA21*rho1(3,3)-GAMMA10*rho1(2,2),-gamma21*rho1(2,3);
                        -gamma20*rho1(3,1), -gamma21*rho1(3,2), -GAMMA21*rho1(3,3)]; %Lindblad term
                    
                    LrhoNEW = [-GAMMA01*rho1(1,1), -gamma01*rho1(1,2), -gamma02*rho1(1,3);
                        -gamma01*rho1(2,1), -GAMMA12*rho1(2,2)+GAMMA01*rho1(1,1),-gamma12*rho1(2,3);                       
                        -gamma02*rho1(3,1), -gamma12*rho1(3,2), GAMMA12*rho1(2,2)] ; %Lindblad term
                    
                    rho = (1-dpol_mat(ind1))*(rho1 + deco*(Lrho+other_rates*LrhoNEW)*dt) + dpol_mat(ind1)*eye(3)/3;   % addition of the rotated density matrix with Lindblad (dissipative) term
                final_time=final_time+dt;
                end
                %           
            end
            final_time_mat(ind1,1)=final_time;
            rhof(:,:,ind3,ind1)=rho;
            temperature_log(ind3,ind1)=Temperature;
        end
    end
   pop0=squeeze(rhof(1,1,:,:));
   pop1=squeeze(rhof(2,2,:,:));
   pop2=squeeze(rhof(3,3,:,:));
   %%%%%%%%%%%%%%
   var_all(:,:,data_ind)=var(pop0);
   sim_plot_data(:,:,data_ind)=pop0;
   
   avg_counter=avg_counter+1;
end
pop0=mean(sim_plot_data,3);
%%
min_sim_all=zeros(n_bombs,avg_number);
max_sim_all=zeros(n_bombs,avg_number);
mean_sim_all=mean(sim_plot_data,1);
var_sim_all=zeros(n_bombs,avg_number);
std_sim_all=zeros(n_bombs,avg_number);
for avg_ind=1:avg_number
    sim_plot_data_ind=sim_plot_data(:,:,avg_ind);  
min_sim_all(:,avg_ind)=squeeze(min(sim_plot_data_ind));
max_sim_all(:,avg_ind)=squeeze(max(sim_plot_data_ind));
var_sim_all(:,avg_ind)=squeeze(var(sim_plot_data_ind));
std_sim_all(:,avg_ind)=squeeze(std(sim_plot_data_ind));
end
% figure
% plot(abs(min_sim_all))
% hold on
% plot(abs(max_sim_all))
% hold on
% plot(squeeze(abs(mean_sim_all)))
% 
% figure
% plot(squeeze(abs(var_sim_all)))
% 
% figure
% plot(squeeze(abs(std_sim_all)))
%%
figure

tiledlayout(2,1)

nexttile
plot(abs(mean(min_sim_all,2)))
hold on
plot(abs(mean(max_sim_all,2)))
hold on
plot(squeeze(abs(mean(mean_sim_all,3))))
title("Simulation: E[p_0] (yellow), Min. or \theta=0 (blue) and Max or theta=pi (red)" )
    xlabel('N');
    ylabel('D_0 outcome');


nexttile
plot(squeeze(abs(mean(std_sim_all(:,1:end),2))))
title("Standard deviation" )
    xlabel('N');
    ylabel('\sigma[p_0]');
%%
% figure
% plot(squeeze(abs(mean(var_sim_all(:,1:end),2))))
%%
% mean_max_sim_all=[0.26;mean(max_sim_all(2:25,1:end),2)];
% figure (103)
% hold on
% plot(abs(mean(min_sim_all,2)))
% hold on
% plot(abs(mean_max_sim_all))
% hold on
% plot(squeeze(abs(mean(mean_sim_all,3))))
%% 
eta_sim=pop0./(pop0+pop2);
pr_sim=pop0./(pop0+pop1);
nr_sim=pop1./(pop0+pop1);


figure
tiledlayout (2,2)

nexttile %1
   Y1=linspace(0,1,n_experiments+1);
   X1=linspace(1,n_bombs,n_bombs);
   [X,Y]=meshgrid(X1,Y1);

   
   surf(X,Y,eta_sim)
   view ([90, -90])
   shading 'flat'
   colorbar
   caxis ([0 1])
   ylim ([0 1])
   xlim ([1 n_bombs])
    title("\eta_c Simulation")
    xlabel('N');
    ylabel('\theta/\pi');
%%%%%
   nexttile %2
   Y1=linspace(0,1,n_experiments+1);
   X1=linspace(1,n_bombs,n_bombs);
   [X,Y]=meshgrid(X1,Y1);
   
   surf(X,Y,abs(pop0))
   view ([90, -90])
   shading 'flat'
   colorbar
   caxis ([0 1])
   ylim ([0 1])
   xlim ([1 n_bombs])

    title("p_0 Simulation")
    xlabel('N');
    ylabel('\theta/\pi');
  %%%%% %3
   nexttile
   Y1=linspace(0,1,n_experiments+1);
   X1=linspace(1,n_bombs,n_bombs);
   [X,Y]=meshgrid(X1,Y1);
   
   surf(X,Y,pr_sim)
   view ([90, -90])
   shading 'flat'
   colorbar
   caxis ([0 1])
   ylim ([0 1])
   xlim ([1 n_bombs])

    title("PR(\theta) Simulation")
    xlabel('N');
    ylabel('\theta/\pi');

    %%%%% %4
   nexttile
   Y1=linspace(0,1,n_experiments+1);
   X1=linspace(1,n_bombs,n_bombs);
   [X,Y]=meshgrid(X1,Y1);
   
   surf(X,Y,nr_sim)
   view ([90, -90])
   shading 'flat'
   colorbar
   caxis ([0 1])
   ylim ([0 1])
   xlim ([1 n_bombs])

    title("NR(\theta) Simulation")
    xlabel('N');
    ylabel('\theta/\pi');
   %%
%    figure
%    plot(final_time_mat)
   %%
    meanpop=mean(pop0);
    varpop=var(pop0);
    popmin=pop0(1,:);
    popmax=pop0(end,:);
if (Random_strengths==0)
 figure
tiledlayout(2,1)
nexttile
 plot(X1,meanpop)
 hold on
 plot(X1,popmin)
 hold on
 plot(X1,popmax)
 title("Simulation: E[p_0] (blue), p_0[\theta=0] (red) and p_0[\theta=\pi] (yellow)" )
    xlabel('N');
    ylabel('D_0 outcome');

 nexttile
  plot(X1,varpop)
  hold on
  plot(squeeze(abs(mean(std_sim_all(:,1:end),2))))
title("Variance and Standard deviation" )
    xlabel('N');
  %  ylabel('\sigma[p_0]');
  %%%%%%%%%%%%%%%%%%%%%%%
elseif (Random_strengths==1)
 figure
tiledlayout(2,1)
nexttile
 plot(X1,meanpop, X1,min(pop0), X1,max(pop0))
 title("Simulation: E[p_0] (blue), Min. (red) and Max. (yellow)" )
    xlabel('N');
    ylabel('D_0 outcome');

 nexttile
  plot(X1,varpop)
  hold on
  plot(squeeze(abs(mean(std_sim_all(:,1:end),2))))
title("Variance and Standard deviation" )
    xlabel('N');
  %  ylabel('\sigma[p_0]');
end

%  figure()
%  hold on
%  plot(X1,varpop)