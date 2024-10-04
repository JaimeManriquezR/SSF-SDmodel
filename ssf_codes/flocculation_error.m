function C = flocculation_error(t, d_error, l_error, c_reference, e_factor, e_components)
arguments
    t
    d_error = 1.0;
    l_error = 1.0;
    c_reference = [0.0026808, 0.0100004, 0, 0.0091, 0.00623, 0.00001998818, 0, 0.00017493]
    e_factor = 2.0;
    e_components = [1, 1, 0, 0, 0, 0, 0, 1];
end

if (d_error <= t) && (t <= d_error + l_error)
    C = (1 + (e_factor/100).*e_components) .* c_reference;
else
    C = c_reference;
end

end