function output = ao_dot(u_mat,v_mat)
    u = squeeze(u_mat);
    v = squeeze(v_mat);
    output = 0;
    for i = 1:length(u)
        output = output + u(i)*v(i);
    end
end
