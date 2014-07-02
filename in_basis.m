function output = in_basis(basis_mat,vec_mat)
    basis = squeeze(basis_mat);
    vec = squeeze(vec_mat);
    inv_basis = inv(basis);
    for i = 1:length(vec)
        output(i) = 0*vec(i);
        for j = 1:length(vec)
            output(i) = output(i) + inv_basis(i,j)*vec(j);
        end
    end
end