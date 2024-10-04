function plot_biofilm(obj, t, options)
arguments
    obj (1,1) SDresults
    t (1,1) = obj.time_final; 
    options.concentrations = sum_concentrations(obj);
    options.text_size = 12;
    options.fig_handle = figure;
end

C_sum = options.concentrations;
text_size = options.text_size;

all_index = (1:length(obj.filter.mesh.cell_centers)).';
water_index = all_index(obj.filter.mesh.cell_centers < 0);
sand_index = [water_index(end); setdiff(all_index, water_index)];

zj_water = obj.filter.mesh.cell_centers(water_index);
zj_sand = obj.filter.mesh.cell_centers(sand_index);

ze_water = zj_water.' + obj.filter.mesh.dz*[-1/2; 1/2];
ze_sand = zj_sand.' + obj.filter.mesh.dz*[-1/2; 1/2];

rhoP = obj.model.particles.density(1);
rhoL = obj.model.liquids.density(1);

phiM = C_sum.("M")/rhoP;
phie = C_sum.("Pe")/rhoP + C_sum.("Le")/rhoL;
%phif = C_sum.("Pf")/rhoP + C_sum.("Lf")/rhoL;

phib = phiM + phie;

t_index = find((obj.frames.time - t) >= 0, 1);

cj_water = phib(water_index, t_index);
cj_sand = phib(sand_index, t_index);

hold on;
plot(kron(cj_water,[1;1]),ze_water(:), ...
    '-','linewidth',1,'markersize',3,'color','#2cde79');
plot(cj_water, zj_water, ...
    'o','linewidth',1,'markersize',3,'color','#2cde79');

plot(kron(cj_sand,[1;1]),ze_sand(:), ...
    '-','linewidth',1,'markersize',3,'color','#fcba03');
plot(cj_sand, zj_sand, ...
    'o','linewidth',1,'markersize',3,'color','#fcba03');

set(gca(),'YDir','reverse');
ylim('tight');
yL = ylabel('Depth $z$ [m]');
xL = xlabel('Biofilm volume fraction $\phi_{\rm b}$');



fL = legend({'$z < 0$: Supernatant water','','$z \geq 0$: Sand bed'});
set([yL, xL, fL], ...
    'interpreter', 'latex', ...
    'FontSize',text_size);
end
