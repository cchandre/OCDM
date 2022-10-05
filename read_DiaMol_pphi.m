function read_DiaMol_pphi
%%
%% Last modified by Cristel Chandre (October 5, 2022)
%% Comments? cristel.chandre@cnrs.fr 
%%

[filename, path] = uigetfile('*.mat');
load([path filename],'data','type_traj')
Nbin = 100;
if strcmp(strtrim(type_traj(1,:)),'all')
    Nd = length(data{1,2}(:,1))/4;
    Nu = length(data{1,3}(:,1))/4;
    Pd = data{1,2}(3*Nd+1:end,end);
    Pu = data{1,3}(3*Nu+1:end,end);
    Pu = Pu(Pu>=150);
    figure, histogram(Pd,Nbin,'Normalization','pdf')
    hold on, histogram(Pu,Nbin,'Normalization','pdf')
else
    N = length(data{1,2}(:,1))/4;
    P = data{1,2}(3*N+1:end,end);
    figure, histogram(P,Nbin,'Normalization','pdf')
end
set(gca,'box','on','FontSize',20,'LineWidth',2)
xlabel('$p_\phi$ (a.u.)','interpreter','latex','FontSize',26)
ylabel('PDF','FontSize',26)

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