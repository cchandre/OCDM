function read_DiaMol_trajectories
%%
%% Last modified by Cristel Chandre (February 21, 2023)
%% Comments? cristel.chandre@cnrs.fr 
%%
choice_repr = 'pphi'; % options are 'pphi' or 'p' (rescaled)
method = 'time_series'; % options are 'phase_diag' or 'time_series' 
variables = [1,4]; % indices of the variables to be displayed (r, phi, pr, pphi)
limit = [];
%close all
[filename, path] = uigetfile('*.mat');
load([path filename],'data','dim','mu','beta');
dim = double(dim);
t = data{1,1};
data = data(~cellfun('isempty',data));
if length(data)==3
    n_t = length(t);
    n_d = length(data{1,2}(:,1))/(2*dim);
    n_b = length(data{1,3}(:,1))/(2*dim);
    y_d = reshape(data{1,2},[n_d,2*dim,n_t]);
    y_b = reshape(data{1,3},[n_b,2*dim,n_t]);
    Pend = y_b(:,4,end);
    y_c = y_b(Pend<=150,:,:);
    n_c = length(y_c(:,1,1));
    y_l = y_b(Pend>=150,:,:);
    n_l = length(y_l(:,1,1));
    if strcmp(method,'time_series')
        if dim==2
            labels = {'$r$ (a.u.)','$\phi$','$p_r$ (a.u.)','$p_\phi$ (a.u.)'};
        else
            labels = {'$r$ (a.u.)','$\theta$','$\phi$','$p_r$ (a.u.)','$p_\theta$ (a.u.)','$p_\phi$ (a.u.)'};
        end
        for it = 1:length(variables)
            subplot(length(variables),1,it)
            plot(t,reshape(y_d(:,variables(it),:),[n_d,n_t]),'r','LineWidth',2)
            hold on
            plot(t,reshape(y_l(:,variables(it),:),[n_l,n_t]),'b','LineWidth',2)
            plot(t,reshape(y_c(:,variables(it),:),[n_c,n_t]),'k','LineWidth',2)
            set(gca,'box','on','FontSize',20,'LineWidth',2)
            xlabel('$t$ (ps)','interpreter','latex','FontSize',26)
            ylabel(labels{variables(it)},'interpreter','latex','FontSize',26)
            xlim([min(t), max(t)])
        end
    elseif strcmp(method,'phase_diag')
        for it = 1:n_t
            if strcmp(choice_repr,'p')
                figure, plot(y_d(:,2,it),(y_d(:,4,it)-mu*y_d(:,1,it).^2*beta*t(it)/2.418884254e-5)./(mu*sqrt(beta)*y_d(:,1,it).^2),'r.')
                hold on
                plot(y_l(:,2,it),(y_l(:,4,it)-mu*y_l(:,1,it).^2*beta*t(it)/2.418884254e-5)./(mu*sqrt(beta)*y_l(:,1,it).^2),'b.')
                plot(y_c(:,2,it),(y_c(:,4,it)-mu*y_c(:,1,it).^2*beta*t(it)/2.418884254e-5)./(mu*sqrt(beta)*y_c(:,1,it).^2),'k.')
                ylabel('$p$','interpreter','latex','FontSize',26)
            elseif strcmp(choice_repr,'pphi')
                figure, plot(y_d(:,2,it),y_d(:,4,it),'r.')
                hold on
                plot(y_l(:,2,it),y_l(:,4,it),'b.')
                plot(y_c(:,2,it),y_c(:,4,it),'k.')
                ylabel('$p_\phi$ (a.u.)','interpreter','latex','FontSize',26)
            end
            hold off
            title(['$t =$ ' num2str(t(it))],'interpreter','latex')
            set(gca,'box','on','FontSize',20,'LineWidth',2)
            xlabel('$\phi$','interpreter','latex','FontSize',26)
            xlim([-pi pi])
            if ~isempty(limit)
                ylim(limit)
            end
            pause(0.5)
            exportgraphics(gcf,'DiaMol.gif','Append',true);
        end
    end
else
    n_t = length(t);
    n = length(data{1,2}(:,1))/(2*dim);
    y = reshape(data{1,2},[n,2*dim,n_t]);
    if strcmp(method,'time_series')
        if dim==2
            labels = {'$r$ (a.u.)','$\phi$','$p_r$ (a.u.)','$p_\phi$ (a.u.)'};
        else
            labels = {'$r$ (a.u.)','$\theta$','$\phi$','$p_r$ (a.u.)','$p_\theta$ (a.u.)','$p_\phi$ (a.u.)'};
        end
        for it = 1:2*dim
            subplot(2*dim,1,it)
            plot(t,reshape(y(:,it,:),[n,n_t]),'b','LineWidth',2)
            set(gca,'box','on','FontSize',20,'LineWidth',2)
            xlabel('$t$ (ps)','interpreter','latex','FontSize',26)
            ylabel(labels{it},'interpreter','latex','FontSize',26)
            xlim([min(t), max(t)])
        end
    elseif strcmp(method,'phase_diag')
        for it = 1:n_t
            if strcmp(choice_repr,'p')
                figure, plot(y(:,2,it),(y(:,4,it)-mu*y(:,1,it).^2*beta*t(it)/2.418884254e-5)./(mu*sqrt(beta)*y(:,1,it).^2),'b.')
                ylabel('$p$','interpreter','latex','FontSize',26)
            elseif strcmp(choice_repr,'pphi')
                figure, plot(y(:,2,it),y(:,4,it),'b.')
                ylabel('$p_\phi$ (a.u.)','interpreter','latex','FontSize',26)
            end
            title(['$t =$ ' num2str(t(it))],'interpreter','latex')
            set(gca,'box','on','FontSize',20,'LineWidth',2)
            xlabel('$\phi$','interpreter','latex','FontSize',26)
            xlim([-pi pi])
            if ~isempty(limit)
                ylim(limit)
            end
            pause(0.5)
            exportgraphics(gcf,'DiaMol.gif','Append',true);
        end
    end
end
%
% Copyright (c) 2022 Cristel Chandre.
% All rights reserved.
%
% Redistribution and use in source and binary forms are permitted provided 
% that the above copyright notice and this paragraph are duplicated in all 
% such forms and that any documentation, advertising materials, and other 
% materials related to such distribution and use acknowledge that the
% software was developed by the CNRS. The name of the CNRS may not be used 
% to endorse or promote products derived from this software without 
% specific prior written permission.

% THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF 
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.