% Evaluate whether Adaptive is a nash equilibria as a function of average
% S_0 and T_0 and Delta.

clf; hold on
m = 10^6;
step = 0.01;
T_range = 0:step:2;
S_range = -1:step:1;
payoff_aa = zeros(length(S_range),length(T_range));
payoff_ac = zeros(length(S_range),length(T_range));
payoff_ad = zeros(length(S_range),length(T_range));
payoff_ca = zeros(length(S_range),length(T_range));
payoff_da = zeros(length(S_range),length(T_range));
aa_ess_games = cell(4,1);
colormap(hot)

a1 = 0;
for delta_payoffs = [0.1, 0.5, 1, 2]
    a1 = a1+1;
    tic
    for i = 1:length(S_range)
        S_0 = S_range(i);

        for j = 1:length(T_range)
            T_0 = T_range(j);

            suckers_payoffs = S_0 - delta_payoffs/2 + delta_payoffs * rand(m,1); %random S
            temptations = T_0 - delta_payoffs/2 + delta_payoffs * rand(m,1); %random T
            
            md = double(suckers_payoffs>0 & temptations<1);
            not_md = 1-md;
            
            qMD = sum(md); %number of md games
            payoff_aa(i,j) = qMD/m; %payoff of adaptive agent vs itself
            payoff_ca(i,j) = qMD/m + (1-qMD/m)*sum(suckers_payoffs.*not_md)/max(m-qMD,1); %average S in games with doc < 0
            payoff_da(i,j) = qMD/m *sum(temptations.*md)/max(qMD,1); %average T in games with doc >= 0
            payoff_ac(i,j) = (1-qMD/m)*sum(temptations.*not_md)/max(m-qMD,1); %average T in games with doc < 0
            payoff_ad(i,j) = qMD/m *sum(suckers_payoffs.*md)/max(qMD,1); %average T in games with doc < 0

        end
    end
    toc

    aa_ess1 = int8(payoff_ca < payoff_aa & payoff_aa > payoff_da);
    aa_ess2 = int8(payoff_ca == payoff_aa & payoff_ac > 1 & payoff_aa > payoff_da);
    aa_ess3 = int8(payoff_da == payoff_aa & payoff_ad > 0 & payoff_ca < payoff_aa);
    aa_ess_games{a1} = aa_ess1;% + aa_ess2 + aa_ess3;

    subplot(2,2,a1);
    surf(T_range,S_range,aa_ess_games{a1},'EdgeColor','none')
    set(gca,'FontSize',12)
    xlabel('\langleT\rangle')
    ylabel('\langleS\rangle')
    title(['\Delta = ',num2str(delta_payoffs)])
    view([0 90]);
    drawnow

end



