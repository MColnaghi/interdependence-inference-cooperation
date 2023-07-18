function doc = DegreeOfCorrespondence(S,T)

actor_control = 0.5 * (1 + S - T);
partner_control = 0.5 * (1 + T - S);
joint_control = 0.5 * (1 - S - T);

doc = (2*actor_control.*partner_control + joint_control.^2) ./ ...
    (actor_control.^2 + partner_control.^2 + joint_control.^2);

end