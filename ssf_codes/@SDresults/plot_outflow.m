function [fig_in, fig_out] = plot_outflow(obj, particles_index, liquids_index, options)
arguments
    obj (1,1) SDresults
    particles_index (1,:) = 1:length(obj.model.particles.Name)
    liquids_index (1,:) = 1:length(obj.model.liquids.Name)
    options.figure_inflow = figure;
    options.figure_outflow = figure;
    options.line_colors = ["r";"g";"b";"m";"k"];
    options.time_interval = [0, Inf];
    options.linestyles = ["-", "-."];
    options.yscale = 'log';
end
time = obj.frames.time;
time(time < options.time_interval(1)) = [];
time(time > options.time_interval(2)) = [];

kP = length(obj.model.particles.Name);

fig_in = options.figure_inflow;
fig_out = options.figure_outflow;
[layout_particle_inflow, layout_liquid_inflow] = deal(tiledlayout(fig_in, length(particles_index) + length(liquids_index), 1));
[layout_particle_outflow, layout_liquid_outflow] = deal(tiledlayout(fig_out, length(particles_index) + length(liquids_index), 1));

in_time = [];
c_inflow = [];

for r = 1:length(obj.simulation_data)
    current_time = obj.simulation_data(r).time(:);
    current_time = current_time(ismember(current_time, time));
    
    c_inflow = [c_inflow; cell2mat(arrayfun(obj.simulation_data(r).inflow_concentrations,setdiff(current_time,in_time),'UniformOutput',false))];

    in_time = unique([in_time; current_time]);
end

%c_inflow = cell2mat(arrayfun(obj.simulation_data.inflow_concentrations,time,'UniformOutput',false));
c_outflow = cell2mat(cellfun(@(z) z(end,:), obj.frames.concentrations{[obj.model.particles.Name; obj.model.liquids.Name],"flowing"}, 'UniformOutput', false)).';
c_outflow = c_outflow(ismember(obj.frames.time, time),:);

particle_inflow = c_inflow(:,particles_index);
particle_outflow = c_outflow(:,particles_index);

liquid_inflow = c_inflow(:,liquids_index + kP);
liquid_outflow = c_outflow(:,liquids_index + kP);

particle_inflow_axes = axes('visible','off').empty;
particle_outflow_axes = axes('visible','off').empty;

liquid_inflow_axes = axes('visible','off').empty;
liquid_outflow_axes = axes('visible','off').empty;

for i = 1:length(particles_index)
    I = particles_index(i);
    particle_inflow_axes(i) = plot_readings(in_time, particle_inflow(:, i), layout_particle_inflow, ...
        DisplayName = obj.model.particles.Name{I}, ...
        Color = options.line_colors(I), ...
        LineStyle = options.linestyles(1));

    particle_outflow_axes(i) = plot_readings(in_time, particle_outflow(:, i), layout_particle_outflow, ...
        DisplayName = obj.model.particles.Name{I}, ...
        Color = options.line_colors(I), ...
        LineStyle = options.linestyles(1));
end

for i = 1:length(liquids_index)
    I = liquids_index(i);
    liquid_inflow_axes(i) = plot_readings(in_time, liquid_inflow(:, i), layout_liquid_inflow, ...
        DisplayName = obj.model.liquids.Name{I}, ...
        Color = options.line_colors(I), ...
        LineStyle = options.linestyles(2));

    liquid_outflow_axes(i) = plot_readings(in_time, liquid_outflow(:, i), layout_liquid_outflow, ...
        DisplayName = obj.model.liquids.Name{I}, ...
        Color = options.line_colors(I), ...
        LineStyle = options.linestyles(2));
end

all_layouts = [layout_particle_inflow, layout_particle_outflow, layout_liquid_inflow, layout_liquid_outflow];

inflow_axes = [particle_inflow_axes, liquid_inflow_axes];
outflow_axes = [particle_outflow_axes, liquid_outflow_axes];
all_axes = [inflow_axes, outflow_axes];

ylim(all_axes,'padded')
xlim(all_axes, 'tight')
yscale(all_axes,options.yscale)

xlabel(all_layouts, "Days")
ylabel(all_layouts, "Concentration [kg/m^3]");
end

function [axes_handle, plot_handle] = plot_readings(t_data, y_data, layout, options)
arguments
    t_data
    y_data
    layout
    options.LineStyle = '-'
    options.LineWidth = 1
    options.DisplayName = ''
    options.Color = 'k'
end
axes_handle = nexttile(layout);
plot_handle = plot(axes_handle, t_data, y_data);
set(plot_handle, ...
    'LineStyle', options.LineStyle, ...
    'LineWidth', options.LineWidth, ...
    'DisplayName', options.DisplayName, ...
    'Color', options.Color);
legend show;
end