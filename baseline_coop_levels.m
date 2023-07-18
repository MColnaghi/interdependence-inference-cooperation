function coop_level = baseline_coop_levels(N,beta_fermi,S_0,T_0)

% Use Markov Chain to calculate probability of fixation of adaptive agents
% inferring the degree of correspondence and cooperating if it is
% non-negative, starting from a population (AA,AllC,AllD) = (N/3,N/3,N/3)

T = baseline_transition_matrix(N,beta_fermi,S_0,T_0);
T_old = T;
T_new = T_old * T;
while sum(sum(abs(T_new - T_old))) > 10^(-9)
    T_old = T_new;
    T_new = T_old * T_old;
end

prob_C_fixates = T_new(:,N+1);

i = N/2 + 1;
coop_level = prob_C_fixates(i);

end