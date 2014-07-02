function norm = ao_norm(vec_mat)
    vec = squeeze(vec_mat);
    norm = 0;
    for i = 1:length(vec)
        norm = norm + vec(i)^2;
    end
    norm = sqrt(norm);
end
