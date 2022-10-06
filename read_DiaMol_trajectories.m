function read_DiaMol_trajectories
%%
%% Last modified by Cristel Chandre (October 6, 2022)
%% Comments? cristel.chandre@cnrs.fr 
%%
close all
[filename, path] = uigetfile('*.mat');
load([path filename],'data','dim','Method','mu','beta');
dim = double(dim);
t = data{1,1};
data = data(~cellfun('isempty',data));
if length(data)==3
    n_t = length(t);
    n_d = length(data{1,2}(:,1))/(2*dim);
    n_b = length(data{1,3}(:,1))/(2*dim);
    y_d = reshape(data{1,2},[n_d,2*dim,n_t]);
    y_b = reshape(data{1,3},[n_b,2*dim,n_t]);
    if strcmp(strtrim(Method),'trajectories')
        if dim==2
            labels = {'$r$','$\phi$','$p_r$','$p_\phi$'};
        else
            labels = {'$r$','$\theta$','$\phi$','$p_r$','$p_\theta$','$p_\phi$'};
        end
        for it = 1:2*dim
            subplot(2*dim,1,it)
            plot(t,reshape(y_d(:,it,:),[n_d,n_t]),'r','LineWidth',2)
            hold on
            plot(t,reshape(y_b(:,it,:),[n_b,n_t]),'b','LineWidth',2)
            set(gca,'box','on','FontSize',20,'LineWidth',2)
            xlabel('$t$ (ps)','interpreter','latex','FontSize',26)
            ylabel(labels{it},'interpreter','latex','FontSize',26)
            xlim([min(t), max(t)])
        end
    elseif strcmp(strtrim(Method),'dissociation')
        for it = 1:n_t
            figure, plot(y_d(:,2,it),y_d(:,4,it)-mu*y_d(:,1:it).^2*beta*t(it),'r.')
            hold on
            figure, plot(y_b(:,2,it),y_b(:,4,it)-mu*y_b(:,1:it).^2*beta*t(it),'b.')
            hold off
            set(gca,'box','on','FontSize',20,'LineWidth',2)
            xlabel('$\phi$','interpreter','latex','FontSize',26)
            ylabel('$\tilde{p}_\phi$',interpreter','latex','FontSize',26)
        end
    end
else
    n_t = length(t);
    n = length(data{1,2}(:,1))/(2*dim);
    y = reshape(data{1,2},[n,2*dim,n_t]);
    if strcmp(strtrim(Method),'trajectories')
        if dim==2
            labels = {'$r$','$\phi$','$p_r$','$p_\phi$'};
        else
            labels = {'$r$','$\theta$','$\phi$','$p_r$','$p_\theta$','$p_\phi$'};
        end
        for it = 1:2*dim
            subplot(2*dim,1,it)
            plot(t,reshape(y(:,it,:),[n,n_t]),'b','LineWidth',2)
            set(gca,'box','on','FontSize',20,'LineWidth',2)
            xlabel('$t$ (ps)','interpreter','latex','FontSize',26)
            ylabel(labels{it},'interpreter','latex','FontSize',26)
            xlim([min(t), max(t)])
        end
    elseif strcmp(strtrim(Method),'dissociation')
        for it = 1:n_t
            figure, plot(y(:,2,it),y(:,4,it)-mu*y(:,1:it).^2*beta*t(it),'b.')
            title(['distribution at $t = $' num2str(t(it)) ' ps'],'interpreter','latex')
            set(gca,'box','on','FontSize',20,'LineWidth',2)
            xlabel('$\phi$','interpreter','latex','FontSize',26)
            ylabel('$\tilde{p}_\phi$','interpreter','latex','FontSize',26)
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