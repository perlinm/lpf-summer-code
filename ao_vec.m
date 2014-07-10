function vector = ao_vec(vec_mat,units)
    vec = squeeze(vec_mat);
    unit_pl = plist('axis','y','yunits',units);
    for i = 1:length(vec)
        vector(i) = ao(vec(i),unit_pl);
    end
end
