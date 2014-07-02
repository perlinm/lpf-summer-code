function vector = ao_vec(vec_mat,units)
    vec = squeeze(vec_mat);
    unit_pl = plist('axis','y','yunits',units);
    vector = [ao(vec(1),unit_pl),ao(vec(2),unit_pl),ao(vec(3),unit_pl)];
end
