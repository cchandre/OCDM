function read_DiaMol_trajectories
%%
%% Last modified by Cristel Chandre (September 13, 2022)
%% Comments? cristel.chandre@cnrs.fr 
%%
close all
[filename, path] = uigetfile('*.mat');
load([path filename]);
dim = double(dim);
t = data{1};
n_t = length(t);
n_d = length(data{2}(:,1))/(2*dim);
n_b = length(data{3}(:,1))/(2*dim);
y_d = reshape(data{2},[n_d,2*dim,n_t]);
y_b = reshape(data{3},[n_b,2*dim,n_t]);
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