function [AA_fixation_probability,coop_level] = markov_chain_exact(N,beta_fermi,matrix,alpha)

% Use Markov Chain to calculate probability of fixation of adaptive agents
% inferring the degree of correspondence and cooperating if it is
% non-negative, starting from a population (AA,AllC,AllD) = (N/3,N/3,N/3)

T = transition_matrix(N,beta_fermi,matrix);
r = numel(T(1,:))-3;

%%% Transform matrix in canonic form
P = [T(:,2:N),T(:,N+2:end-1),T(:,1),T(:,N+1),T(:,end)]; % rearrange rows
P = [P(2:N,:); P(N+2:end-1,:); P(1,:); P(N+1,:); P(end,:)]; %  rearrange columns

Q = P(1:end-3,1:end-3);
R = P(1:end-3,end-2:end);
%N_inv = inv(eye(r)-Q); 
%B1 = N_inv*R;
B = (eye(r)-Q)\R;

prob_A_fixates = B(:,end);
prob_C_fixates = B(:,2);

% prob_A_fixates_old = T_new(:,end);
% prob_C_fixates_old = T_new(:,N+1);

% prob_A_fixates_old = [T_new(2:N,end);T_new(N+2:end-1,end)];
% prob_C_fixates_old = [T_new(2:N,N+1);T_new(N+2:end-1,N+1)];
%prob_D_fixates = T_new(:,1);


X = round(N/3);
Y = round(N/3);
num_list = N+1 - (0:(X-1));
i = sum(num_list) + Y + 1;

% disp([prob_A_fixates(i-2),prob_A_fixates_old(i),prob_C_fixates(i-2),prob_C_fixates_old(i)])


AA_fixation_probability = prob_A_fixates(i);
coop_level = AA_fixation_probability*alpha + prob_C_fixates(i);

end