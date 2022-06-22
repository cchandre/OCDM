function read_DiaMol_dissociation
%%
%% Last modified by Cristel Chandre (June 20, 2022)
%% Comments? cristel.chandre@cnrs.fr 
%%

[filename, path] = uigetfile('*.txt');
data = importdata([path filename],' ',3) ;
[t, index] = sort(data.data(:,1));
proba = data.data(index,2);
figure, plot(t,proba,'LineWidth',3)
set(gca,'box','on','FontSize',20,'LineWidth',2)
xlabel('$t$ (ps)','interpreter','latex','FontSize',26)
ylabel('probability','FontSize',26)

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