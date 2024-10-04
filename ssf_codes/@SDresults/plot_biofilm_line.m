function bf_line = plot_biofilm_line(obj, t, options)
arguments
    obj (1,1) SDresults
    t (1,1) = obj.time_final; 
    options.concentrations = sum_concentrations(obj);
    options.text_size = 12;
    options.fig_handle = figure;
    options.color = 'k';
    options.linestyle = '-';
    options.linewidth = 1;
    options.display_name = "biofilm";
end

C_sum = options.concentrations;
text_size = options.text_size;

rhoP = obj.model.particles.density(1);
rhoL = obj.model.liquids.density(1);

phiM = C_sum.("M")/rhoP;
phie = C_sum.("Pe")/rhoP + C_sum.("Le")/rhoL;
%phif = C_sum.("Pf")/rhoP + C_sum.("Lf")/rhoL;

phib = phiM + phie;

bf_axes = axes(options.fig_handle, 'NextPlot', 'add');

t_index = find((obj.frames.time - t) >= 0, 1);

bf_line = plot(bf_axes, phib(:, t_index), obj.filter.mesh.cell_centers, ...
     options.linestyle, ...
    'linewidth',options.linewidth, ...
    'color',options.color, ...
    'DisplayName', options.display_name, ...
    'MarkerIndices',1:100:1001);
grid(bf_axes, "on")

set(bf_axes,'YDir','reverse');
ylim(bf_axes, obj.filter.domain);
yL = ylabel(bf_axes, 'Depth $z$ [m]');
xL = xlabel(bf_axes, 'Biofilm volume fraction $\phi_{\rm b}$');
fL = legend(bf_axes, 'show');
set([yL, xL, fL], ...
    'interpreter', 'latex', ...
    'FontSize',text_size);
end
