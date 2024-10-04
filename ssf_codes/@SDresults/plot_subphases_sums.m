function plot_subphases_sums(obj, options)
arguments
    obj (1,1) SDresults
    options.sub_phases_index = {'M'; 'Pe'; 'Le'; 'Pf'; 'Lf'};
    options.fig_handle = figure;
end
sub_phases = options.sub_phases_index;

phase_layout = tiledlayout(options.fig_handle, 'flow');
for s = 1:length(sub_phases)
    %figure
    p_axes(s) = nexttile;
    mesh(p_axes(s),T,Z,C_sum.(sub_phases{s}))
    title(p_axes(s), sprintf("$\\displaystyle \\sum_{i = 1}^{k_{\\rm %s}} C^i_{\\rm %s}$",sub_phases{s},sub_phases{s}),'Interpreter','latex')
    %colormap summer
end
set(p_axes,'XLimitMethod','tight')
set(p_axes,'YLim',obj.filter.domain)
set(p_axes, 'YDir', 'reverse')
view(p_axes, 2)
title(phase_layout, "(Sum of) Global concentrations by sub-phase")
end
