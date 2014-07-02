function output = ao_cross(u_mat,v_mat)
    u = squeeze(u_mat);
    v = squeeze(v_mat);
    output = [u(2).*v(3) - u(3).*v(2),...
              u(3).*v(1) - u(1).*v(3),...
              u(1).*v(2) - u(2).*v(1)];
end
