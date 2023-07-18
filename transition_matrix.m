function T = transition_matrix(N,beta_fermi,m)

a11 = m(1,1);
a12 = m(1,2);
a13 = m(1,3);
a21 = m(2,1);
a22 = m(2,2);
a23 = m(2,3);
a31 = m(3,1);
a32 = m(3,2);
a33 = m(3,3);

num_list = N+1 - (0:(N-1));
T = zeros(sum(num_list)+1);
%T = sparse([]);
sum_k = 0;

T(1,1) = 1;
X = 0;

for Y = 1:N
    Z =(N-X-Y);
    payoff_A = (a11*(X-1) + a12*Y + a13*Z)/(N-1);
    payoff_C = (a21*X + a22*(Y-1) + a23*Z)/(N-1);
    payoff_D = (a31*X + a32*Y + a33*(Z-1))/(N-1);

    i = Y + 1;
    T(i,i - 1) = (N-X-Y)/N * Y/(N-1)/ (1 + exp(-beta_fermi*(payoff_D - payoff_C)));
    T(i,i + 1) = Y/N * (N-X-Y)/(N-1)/ (1 + exp(-beta_fermi*(payoff_C - payoff_D)));
    T(i,i + N-X+1) = X/N * (N-X-Y)/(N-1)/ (1 + exp(-beta_fermi*(payoff_A - payoff_D)));
    T(i,i + N-X) = X/N * Y/(N-1)/ (1 + exp(-beta_fermi*(payoff_A - payoff_C)));
    T(i,i) = 1 - T(i,i - 1) -T(i,i + 1) - T(i,i + N-X+1) - T(i,i + N-X);
end
sum_k = sum_k + (N-X+1);

for X = 1:N-1
    Y = 0;
    Z =(N-X-Y);
    payoff_A = (a11*(X-1) + a12*Y + a13*Z)/(N-1);
    payoff_C = (a21*X + a22*(Y-1) + a23*Z)/(N-1);
    payoff_D = (a31*X + a32*Y + a33*(Z-1))/(N-1);

    i = sum_k + 1;
    T(i,i - N-2+X) = (N-X-Y)/N * X/(N-1)/ (1 + exp(-beta_fermi*(payoff_D - payoff_A)));
    T(i,i - N-1+X) = Y/N * X/(N-1)/ (1 + exp(-beta_fermi*(payoff_C - payoff_A)));
    T(i,i + 1) = Y/N * (N-X-Y)/(N-1)/ (1 + exp(-beta_fermi*(payoff_C - payoff_D)));
    T(i,i + N-X+1) = X/N * (N-X-Y)/(N-1)/ (1 + exp(-beta_fermi*(payoff_A - payoff_D)));
    T(i,i) = 1 -T(i,i - N-1+X) -T(i,i - N-2+X) -T(i,i + 1) - T(i,i + N-X+1);
    
    for Y = 1:(N-X-1)
    Z =(N-X-Y);
    payoff_A = (a11*(X-1) + a12*Y + a13*Z)/(N-1);
    payoff_C = (a21*X + a22*(Y-1) + a23*Z)/(N-1);
    payoff_D = (a31*X + a32*Y + a33*(Z-1))/(N-1);
        
        i = sum_k + Y + 1;
    T(i,i - N-2+X) = (N-X-Y)/N * X/(N-1)/ (1 + exp(-beta_fermi*(payoff_D - payoff_A)));
    T(i,i - N-1+X) = Y/N * X/(N-1)/ (1 + exp(-beta_fermi*(payoff_C - payoff_A)));
    T(i,i - 1) = (N-X-Y)/N * Y/(N-1)/ (1 + exp(-beta_fermi*(payoff_D - payoff_C)));
    T(i,i + 1) = Y/N * (N-X-Y)/(N-1)/ (1 + exp(-beta_fermi*(payoff_C - payoff_D)));
    T(i,i + N-X+1) = X/N * (N-X-Y)/(N-1)/ (1 + exp(-beta_fermi*(payoff_A - payoff_D)));
    T(i,i + N-X) = X/N * Y/(N-1)/ (1 + exp(-beta_fermi*(payoff_A - payoff_C)));
    T(i,i) = 1 -T(i,i - N-2+X) -T(i,i - N-1+X) - T(i,i - 1)...
                -T(i,i + 1) - T(i,i + N-X+1) - T(i,i + N-X);
    end
    
    Y = N-X;
    Z =(N-X-Y);
    payoff_A = (a11*(X-1) + a12*Y + a13*Z)/(N-1);
    payoff_C = (a21*X + a22*(Y-1) + a23*Z)/(N-1);
    
    i = sum_k + Y + 1;
    T(i,i - N-2+X) = (N-X-Y)/N * X/(N-1)/ (1 + exp(-beta_fermi*(payoff_D - payoff_A)));
    T(i,i - N-1+X) = Y/N * X/(N-1)/ (1 + exp(-beta_fermi*(payoff_C - payoff_A)));
    T(i,i - 1) = (N-X-Y)/N * Y/(N-1)/ (1 + exp(-beta_fermi*(payoff_D - payoff_C)));
    T(i,i + N-X) = X/N * Y/(N-1)/ (1 + exp(-beta_fermi*(payoff_A - payoff_C)));
    T(i,i) = 1 -T(i,i - N-2+X) -T(i,i - N-1+X) - T(i,i - 1)...
                - T(i,i + N-X);

    sum_k = sum_k + (N-X+1);
end

T(end,end)=1;
end







