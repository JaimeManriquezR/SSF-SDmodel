function [f1, f2, f3] = plot_biofilm_evolution(obj, time_vector, options)
arguments
    obj (1,1) SDresults
    time_vector (:,1) = obj.time_final; 
    options.concentrations = sum_concentrations(obj);
    options.text_size = 12;
    options.fig_full_handle = figure;
    options.fig_zoom_handle = figure;
    options.fig_change_handle = figure;
    options.zlim = 'tight';
end

marker_list = string(('-.osvx.-v+*.xsdv><').');
color_list = {'k', 'm', 'r', 'b', "#77AC30", "#ffa500"};

f1 = options.fig_full_handle;
f2 = options.fig_zoom_handle;
f3 = options.fig_change_handle;

C_sum = options.concentrations;
text_size = options.text_size;

zj  = obj.filter.mesh.cell_centers;
ze = zj.' + obj.filter.mesh.dz*[-1/2; 1/2];

rhoP = obj.model.particles.density(1);
rhoL = obj.model.liquids.density(1);

phiM = C_sum.("M")/rhoP;
phie = C_sum.("Pe")/rhoP + C_sum.("Le")/rhoL;
%phif = C_sum.("Pf")/rhoP + C_sum.("Lf")/rhoL;

phib = phiM + phie;

bf_axes = axes(f1, 'NextPlot', 'add');


for i = 1:length(time_vector)
t = time_vector(i);

t_index = find((obj.frames.time - t) >= 0, 1);

cj = phib(:, t_index);

% plot(bf_axes, kron(cj,[1;1]),ze(:), ...
%     '-','linewidth',1,'markersize',3,'color',color_list{i}, ...
%     'DisplayName', sprintf("Day $%i$", round(t)));

plot(bf_axes, cj, zj(:), ...
    '-' + marker_list(i),'linewidth',1,'markersize',6,'color',color_list{i}, ...
    'DisplayName', sprintf("Day $%i$", round(t)), 'MarkerIndices', 1:100:1001);
end
grid(bf_axes, "on")

set(bf_axes,'YDir','reverse');
ylim(bf_axes, obj.filter.domain);
yL = ylabel(bf_axes, 'Depth $z$ [m]');
xL = xlabel(bf_axes, 'Biofilm volume fraction $\phi_{\rm b}$');
fL = legend(bf_axes, 'show');
set([yL, xL, fL], ...
    'interpreter', 'latex', ...
    'FontSize',text_size);

zoom_axes = copyobj(f1.Children, f2);
ylim(zoom_axes(2), options.zlim)
set(f2.Children(2).Children,'MarkerIndices',1:10:1001)

cj0 = phib(:, 1);
for i = 1:floor(obj.time_final)

t_index = find((obj.frames.time - i) >= 0, 1);

cj1 = phib(:, t_index);

dcdt(i) = 100 * ((norm(cj1,2) - norm(cj0,2))/norm(cj1,2));
cj0 = cj1;
end

dcdt_axes = axes(f3);
plot(dcdt_axes, 1:floor(obj.time_final), dcdt, 'ro--')

grid(dcdt_axes,"on")

yyL = ylabel(dcdt_axes, 'Relative change [\%]');
xxL = xlabel(dcdt_axes, 'Days');
yscale(dcdt_axes,'log')
set([yyL, xxL], ...
    'interpreter', 'latex', ...
    'FontSize',text_size);
ylim(dcdt_axes, 'tight');
xlim(dcdt_axes, 'tight')

end