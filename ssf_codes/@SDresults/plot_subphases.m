function plot_subphases(obj, options)
arguments
    obj (1,1) SDresults
    options.sub_phases = [{["particles","matrix"];
            ["particles", "enclosed"];
            ["liquids", "enclosed"];
            ["particles", "flowing"];
            ["liquids", "flowing"]}, {'M'; 'Pe'; 'Le'; 'Pf'; 'Lf'}];
    %options.fig_handle = figure;
end
sub_phases = options.sub_phases;

for s = 1:length(sub_phases(:,2))
    sp_axes = [];
    component_type = sub_phases{s,1}(1);
    volume_type = sub_phases{s,1}(2);
    phase_name = sub_phases{s,2};

    figure
    component_layout(s) = tiledlayout('flow');
    title(component_layout(s), "Global concentrations")
    subtitle(component_layout(s),...
        sprintf("$\\mbox{\\boldmath $C$}_{\\rm %s}$",phase_name), ...
        'Interpreter', 'latex');
    for i = 1:length(obj.model.(component_type).Name)
        Ci_gamma = obj.frames.concentrations{obj.model.(component_type).Name{i},volume_type}{:};
        sp_axes(i) = nexttile;
        mesh(sp_axes(i), T,Z,Ci_gamma);
        title(sp_axes(i), sprintf("$\\rm %s$",obj.model.(component_type).Name{i}), ...
            'Interpreter','latex')
    end
    set(sp_axes,'XLimitMethod','tight')
    set(sp_axes,'YLim',obj.filter.domain)
    set(sp_axes, 'YDir', 'reverse')
    view(sp_axes, 2)
end
end