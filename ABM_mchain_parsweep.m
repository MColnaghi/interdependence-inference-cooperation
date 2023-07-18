% Uses markov chain to calculate the fixation probability of adaptive
% agents for various interdependence distributions with average game
% g = (1, S_0, T_0, 0).

clc
%figure; hold on
clf; hold on

c = load('colMap.mat');
colMap = c.colMap([1:2:80,81:176,176:2:end],:)/255;
%mex Shuffle.c

rng(1,'simdTwister');

%% Simulation parameters

N = 60; % Population size
m = 10^6; % Number of randomly sampled games
delta_payoffs = 2; % Width of s and t range
beta_fermi = 1; %'Temperature' of selection for Moran process; 0 = random drift; inf = hard selection
cost = 0.1; % cost of inference

S_0_range = -1:0.02:1;
T_0_range = 0:0.02:2;

AA_fixation_prob1 = zeros(length(S_0_range),length(T_0_range));
coop_levels1 = zeros(length(S_0_range),length(T_0_range));
AA_fixation_prob2 = zeros(length(S_0_range),length(T_0_range));
coop_levels2 = zeros(length(S_0_range),length(T_0_range));
baseline_coop = zeros(length(S_0_range),length(T_0_range));

deg_of_interdep = zeros(numel(AA_fixation_prob1),1);
fixprob_1 = zeros(numel(AA_fixation_prob1),1);
fixprob_2 = zeros(numel(AA_fixation_prob1),1);
k1=0;

%% Calculations
tic
for i = 1:numel(S_0_range)
    S_0 = S_0_range(i);
    disp(['Round ',num2str(1+(i-1)*numel(T_0_range)),' of ',...
        num2str(numel(S_0_range)*numel(T_0_range)),'. Time Elapsed = ',...
        num2str(round(toc/3600,2)), ' hours'])
    
    for j = 1:numel(T_0_range)
        k1=k1+1;
        T_0 = T_0_range(j);
        suckers_payoffs = S_0- delta_payoffs/2 + delta_payoffs * rand(m,1); % generate random s.payoffs
        temptations = T_0- delta_payoffs/2 + delta_payoffs * rand(m,1); % generate random t. to defect
        
        baseline_coop(i,j) = baseline_coop_levels(N,beta_fermi,S_0,T_0);
        
        % DOC-classification adaptive agents
        [A1,alpha] = adaptive_payoffs_doc_class(cost,S_0,T_0,suckers_payoffs,temptations);
        [AA_fixation_prob1(i,j), coop_levels1(i,j)] = markov_chain_exact(N,beta_fermi,A1,alpha);
        
        % DOC-probability adaptive agents
        [A2,gamma] = adaptive_payoffs_game_class(cost,S_0,T_0,suckers_payoffs,temptations);
        [AA_fixation_prob2(i,j), coop_levels2(i,j)] = markov_chain_exact(N,beta_fermi,A2,gamma);
        
        deg_of_interdep(k1) = degree_of_interdependence([S_0,T_0]);
        fixprob_1(k1) = AA_fixation_prob1(i,j);
        fixprob_2(k1) = AA_fixation_prob2(i,j);
        
    end
end
toc

%% Plot
clf; hold on
ax1 = subplot(2,2,1); hold on; title('A')
ax2 = subplot(2,2,2); hold on; title('B')
ax3 = subplot(2,2,3); hold on; title('C')
ax4 = subplot(2,2,4); hold on; title('D')
colormap(ax1,hot)
colormap(ax2,hot)
% colormap(ax3,turbo)
% colormap(ax4,turbo)
colormap(ax3,jet)
colormap(ax4,jet)

surf(ax1,T_0_range,S_0_range,AA_fixation_prob1,'EdgeColor','none')
surf(ax2,T_0_range,S_0_range,AA_fixation_prob2,'EdgeColor','none')
surf(ax3,T_0_range,S_0_range,coop_levels1 -baseline_coop,'EdgeColor','none')
surf(ax4,T_0_range,S_0_range,coop_levels2 -baseline_coop,'EdgeColor','none')

clim(ax1,[0,1])
clim(ax2,[0,1])
clim(ax3,[-0.3,.6])
clim(ax4,[-0.3,.6])

for ax = [ax1,ax2,ax3,ax4]
    colorbar(ax)
    xlabel(ax,'\langleT\rangle')
    ylabel(ax,'\langleS\rangle')
    set(ax,'Fontsize',12)
    view(ax,[0 90]);
    plot3(ax,[1,1],[-1,1],[1,1],'linewidth',1,'color','k')
    plot3(ax,[0,2],[0,0],[1,1],'linewidth',1,'color','k')
end

%%
wid = 0.31;
hei = 0.35;
ax1.Position = [0.09    0.6    wid    hei];
ax2.Position = [0.59    0.6    wid    hei];
ax3.Position = [0.09    0.11    wid    hei];
ax4.Position = [0.59    0.11    wid    hei];
