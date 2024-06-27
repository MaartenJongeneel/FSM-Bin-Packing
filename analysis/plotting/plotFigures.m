set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex'); set(groot,'defaultLegendInterpreter','latex');
close all;

color.sim = [201, 174, 145]/255;
color.sim_hard = [220,217,208]/255;
color.exp = [51,92,103]/255;

doSave = true;

%Image size
xsize = 400;
ysize = 350;

%%
figure('rend','painters','Position',[400,260,xsize,ysize]);
ha = tight_subplot(4,2,[.05 .1],[.08 .013],[0.1 0.01]);  %[gap_h gap_w] [lower upper] [left right]
for ii=0:2
    axes(ha(ii*2+1));
    plot(t_vec,pos_exp(ii+1,:),'.-','color',color.exp,'linewidth',1,'MarkerSize',4); hold on;
%     plot(t_vec,pos_sim_hard(ii+1,indx_hf2lf(indx_opt,:)),'.-','color',color.sim_hard);
    plot(t_vec,pos_sim(ii+1,:),'.-','color',color.sim,'linewidth',1,'MarkerSize',4);
    % plot(t_vec,pose_exp_nf(ii+1,:),'.-','color','g');
    if length(vlines) > 1
        for jj=1:size(vlines,2)
            xline(vlines(1,jj),'k-.');
        end
    end    
    grid on
%     xlabel('time [s]')
    ylabel(['$p_' ax(ii+1) '$ [m]'])
    xlim(xlims)
    

    axes(ha((ii+1)*2));
    plot(t_vec,rad2deg(rot_exp(ii+1,:)),'.-','color',color.exp,'linewidth',1,'MarkerSize',4); hold on;
%     plot(t_vec,rad2deg(rot_sim_hard(ii+1,indx_hf2lf(indx_opt,:))),'.-','color',color.sim_hard);
    plot(t_vec,rad2deg(rot_sim(ii+1,:)),'.-','color',color.sim,'linewidth',1,'MarkerSize',4);
    % plot(t_vec,rad2deg(rot_exp_nf(ii+1,:)),'.-','color','g');
    if length(vlines) > 1
        for jj=1:size(vlines,2)
            xline(vlines(1,jj),'k-.');
        end
    end
    grid on
%     xlabel('time [s]')
    ylabel(['$\theta_' ax(ii+1) '$ [deg]'])
    xlim(xlims)
end

%     subplot(4,2,7)
    axes(ha(7));
%     plot(t_vec,err_pos_hard,'.-','color',color.sim_hard); hold on
    plot(t_vec,err_pos,'.-','color',color.sim,'linewidth',1,'MarkerSize',4); 
    hold on
    grid on
    if length(vlines) > 1
        for jj=1:size(vlines,2)
            xline(vlines(1,jj),'k-.');
        end
    end
    xlabel('time [s]')
    ylabel('pos\_error [m]')
    xlim(xlims)
    yticks([0 0.025 0.05])

%     subplot(4,2,8)
    axes(ha(8)); 
%     plot(t_vec,rad2deg(err_rot_hard),'.-','color',color.sim_hard); hold on;
    plot(t_vec,rad2deg(err_rot),'.-','color',color.sim,'linewidth',1,'MarkerSize',4);
    hold on
    grid on
    if length(vlines) > 1
        for jj=1:size(vlines,2)
            xline(vlines(1,jj),'k-.');
        end
    end
    xlabel('time [s]')
    ylabel('rot\_error [deg]')
    xlim(xlims)

    if doSave ==1; fig = gcf; fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,append('Seq1_Box016_poses.pdf'),'-dpdf','-painters')
    end

%%
figure('rend','painters','Position',[900,260,xsize,ysize]);
ha = tight_subplot(4,2,[.05 .1],[.08 .013],[0.1 0.01]);  %[gap_h gap_w] [lower upper] [left right]
axes(ha(1));
% subplot(4,2,1)
plot(t_real,ft_real(1,:),'.-','color',color.exp,'linewidth',1,'MarkerSize',4);hold on;
% plot(t_sim,ft_sim_hard(1,hf2lf(indx_opt_ft,:)),'.-','color',color.sim_hard);
plot(t_sim,ft_sim(1,:),'.-','color',color.sim,'linewidth',1,'MarkerSize',4);grid on;
if length(vlines) > 1
    for jj=1:size(vlines,2)
        xline(vlines(1,jj),'k-.');
    end
end 
ylabel('$f_x [N]$')
% xlabel('Time [s]')
xlim(xlims)
ylim([-10 12])

axes(ha(3));
% subplot(4,2,3)
plot(t_real, ft_real(2,:),'.-','color',color.exp,'linewidth',1,'MarkerSize',4);hold on;
% plot(t_sim,ft_sim_hard(2,hf2lf(indx_opt_ft,:)),'.-','color',color.sim_hard);
plot(t_sim, ft_sim(2,:),'.-','color',color.sim,'linewidth',1,'MarkerSize',4);grid on;

if length(vlines) > 1
    for jj=1:size(vlines,2)
        xline(vlines(1,jj),'k-.');
    end
end
ylabel('$f_y [N]$')
% xlabel('Time [s]')
xlim(xlims)
ylim([-20 20])

axes(ha(5));
% subplot(4,2,5)
plot(t_real,ft_real(3,:),'.-','color',color.exp,'linewidth',1,'MarkerSize',4);hold on;
% plot(t_sim,ft_sim_hard(3,hf2lf(indx_opt_ft,:)),'.-','color',color.sim_hard);
plot(t_sim,ft_sim(3,:),'.-','color',color.sim,'linewidth',1,'MarkerSize',4);grid on;

if length(vlines) > 1
    for jj=1:size(vlines,2)
        xline(vlines(1,jj),'k-.');
    end
end
ylabel('$f_z [N]$')
% xlabel('Time [s]')
xlim(xlims)
ylim([-50 50])

axes(ha(7));
% subplot(4,2,7)
% plot(t_sim,f_err_hard,'.-','color',color.sim_hard); hold on;
plot(t,f_err,'.-','color',color.sim,'linewidth',1,'MarkerSize',4);grid on; 
if length(vlines) > 1
    for jj=1:size(vlines,2)
        xline(vlines(1,jj),'k-.');
    end
end
ylabel('$e_f [N]$')
xlabel('Time [s]')
xlim(xlims)

axes(ha(2));
% subplot(4,2,2)
plot(t_real,ft_real(4,:),'.-','color',color.exp,'linewidth',1,'MarkerSize',4);hold on;
% plot(t_sim,ft_sim_hard(4,hf2lf(indx_opt_ft,:)),'.-','color',color.sim_hard); hold on;
plot(t_sim,ft_sim(4,:),'.-','color',color.sim,'linewidth',1,'MarkerSize',4);
if length(vlines) > 1
    for jj=1:size(vlines,2)
        xline(vlines(1,jj),'k-.');
    end
end
ylabel('$\tau_x [Nm]$')
% xlabel('Time [s]')
xlim(xlims)
ylim([-5 5])

axes(ha(4));
% subplot(4,2,4)
plot(t_real, ft_real(5,:),'.-','color',color.exp,'linewidth',1,'MarkerSize',4);hold on;
% plot(t_sim,ft_sim_hard(5,hf2lf(indx_opt_ft,:)),'.-','color',color.sim_hard); hold on;
plot(t_sim, ft_sim(5,:),'.-','color',color.sim,'linewidth',1,'MarkerSize',4);grid on;

if length(vlines) > 1
    for jj=1:size(vlines,2)
        xline(vlines(1,jj),'k-.');
    end
end
ylabel('$\tau_y [Nm]$')
% xlabel('Time [s]')
xlim(xlims)
ylim([-2.5 4.5])

axes(ha(6));
% subplot(4,2,6)
plot(t_real,ft_real(6,:),'.-','color',color.exp,'linewidth',1,'MarkerSize',4);hold on;
% plot(t_sim,ft_sim_hard(6,hf2lf(indx_opt_ft,:)),'.-','color',color.sim_hard);
plot(t_sim,ft_sim(6,:),'.-','color',color.sim,'linewidth',1,'MarkerSize',4);grid on;
if length(vlines) > 1
    for jj=1:size(vlines,2)
        xline(vlines(1,jj),'k-.');
    end
end
ylabel('$\tau_z [Nm]$')
% xlabel('Time [s]')
xlim(xlims)
ylim([-0.2 2])

axes(ha(8));
% subplot(4,2,8)
% plot(t_sim,t_err_hard,'.-','color',color.sim_hard); hold on;
plot(t,t_err,'.-','color',color.sim,'linewidth',1,'MarkerSize',4);grid on; hold on;
if length(vlines) > 1
    for jj=1:size(vlines,2)
        xline(vlines(1,jj),'k-.');
    end
end
ylabel('$e_{\tau} [Nm]$')
xlabel('Time [s]')
xlim(xlims)

    if doSave ==1; fig = gcf; fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,append('Seq1_Box016_wrenches.pdf'),'-dpdf','-painters')
    end

% %%
% 
% %% Wrenches
% 
% t_loop = -0.02:0.001:0.00;
% 
% [hf2lf,indx_opt_ft]=TimeAlign(t_loop,t_sim,t_sim_hard,ft_sim(1,:)',ft_sim_hard(1,:)',125);
% 
% f_err_hard = vecnorm(ft_real(1:3,:)-ft_sim_hard(1:3,hf2lf(indx_opt_ft,:)));
% t_err_hard = vecnorm(ft_real(4:6,:)-ft_sim_hard(4:6,hf2lf(indx_opt_ft,:)));
% 
% 
% %% Poses
% t_loop = -0.05:0.001:0.05;
% 
% [indx_hf2lf,indx_opt]=TimeAlign(t_loop,t_vec,t_sim_hard,pos_sim(1,:)',pos_sim_hard(1,:)',125);
% 
% err_pos_hard = vecnorm(pos_exp-pos_sim_hard(:,indx_hf2lf(indx_opt,:)));
% err_rot_hard = vecnorm(rot_exp-rot_sim_hard(:,indx_hf2lf(indx_opt,:)));
% %%
% function [indx_hf2lf,indx_opt]=TimeAlign(t_vec,t_lf,t_hf,data_lf,data_hf,Fhf)
% %t_vec : Time vector over which we vary
% %t_lf  : Low-Frequency time signal (nx1)
% %t_hf  : High-Frequency time signal (mx1)
% %data_lf : Low Frequency data signal (nx1)
% %data_hf : High Frequency data signal (mx1)
% %Fhf     : Frequency of High Frequency signal (1x1)
% %
% %indx_hf2lf : index showing which hf indexes are used for each lf index (for each element of t_vec)
% %indx_opt : index of t_vec that gives the optimum value
% cnt=1;
% for ii = 1:length(t_vec)
%     tel=1;
%     for i = 1:length(t_lf)
%         abs_dif = abs(t_lf(i) - (t_hf + t_vec(ii))); %Compute the time difference
%         [~,indx_hf2lf(ii,tel)] = min(abs_dif); %Find the indices for which the time difference is minimal
%         tel=tel+1;
%         cnt=cnt+1;
%         textwaitbar(cnt, length(t_lf)*length(t_vec),'Computing time alignment')
%     end
%     indx_con(:,ii) = abs(t_lf-(t_hf(indx_hf2lf(ii,:))+t_vec(ii))) < (1/Fhf)/2; %Do not consider indices for which the time difference is more than half the sampling time
%     Nerror(ii) = norm(data_lf(indx_con(:,ii)) - data_hf(indx_hf2lf(ii,indx_con(:,ii)))); %For the indices we do consider, compute the error
% end
% [~,indx_opt] = min(Nerror); %Find the time index for which the error is minimal
% end