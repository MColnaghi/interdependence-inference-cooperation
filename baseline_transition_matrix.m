function T = baseline_transition_matrix(N,beta_fermi,S_0,T_0)


T = zeros(N+1);

T(1,1) = 1;

for i = 2:N
    X = i-1;
    Y = N-X;
    payoff_C = ((X-1) + S_0*Y)/(N-1);
    payoff_D = T_0*X/(N-1);
    T(i,i+1) = (X/N)*Y/(N-1)/(1 + exp(-beta_fermi*(payoff_C-payoff_D)));
    T(i,i-1) = X/(N-1)*Y/N/(1 + exp(-beta_fermi*(payoff_D-payoff_C)));
    T(i,i) = 1 - T(i,i+1) - T(i,i-1);
end

T(end,end) = 1;

end

