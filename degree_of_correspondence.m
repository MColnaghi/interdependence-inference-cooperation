function deg_of_cor = degree_of_correspondence(G)
actor_control = 0.5*(1 + G(1) - G(2));
partner_control = 0.5*(1 - G(1) + G(2));
joint_control = 0.5*(1 - G(1) - G(2));
deg_of_cor = (2*actor_control.*partner_control + joint_control.^2)./...
    (actor_control.^2 + partner_control.^2 + joint_control.^2);
end