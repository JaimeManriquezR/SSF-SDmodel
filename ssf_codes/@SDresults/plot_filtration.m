function plot_filtration(obj, particles_index, liquids_index, options)
arguments
    obj (1,1) SDresults
    particles_index (1,:) = 1:length(obj.model.particles.Name)
    liquids_index (1,:) = 1:length(obj.model.liquids.Name)
    options.fig_handle = figure;
    options.line_colors = ["r";"g";"b";"m";"k"];
    options.time_interval = [0, Inf];
end  
time = obj.frames.time;
time(time < options.time_interval(1)) = [];
time(time > options.time_interval(2)) = [];

kP = length(obj.model.particles.Name);

fig_handle = options.fig_handle;
[layout_particle, layout_liquid] = deal(tiledlayout(fig_handle, length(particles_index) + length(liquids_index), 1));
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
particle_removal = 1 - (particle_outflow ./ particle_inflow);

liquid_inflow = c_inflow(:,liquids_index + kP);
liquid_outflow = c_outflow(:,liquids_index + kP);
liquid_removal = 1 - (liquid_outflow ./ liquid_inflow);

for i = 1:length(particles_index)
    particle_removal_axes(i) = plot_readings(in_time, 100 * particle_removal(:, i), layout_particle, ...
        DisplayName = obj.model.particles.Name{i}, ...
        Color = options.line_colors(i));
end

for i = 1:length(liquids_index)
    liquid_removal_axes(i) = plot_readings(in_time, 100 * liquid_removal(:, i), layout_liquid, ...
        DisplayName = obj.model.liquids.Name{i}, ...
        Color = options.line_colors(i));
end

all_layouts = [layout_particle, layout_liquid];

all_axes = [particle_removal_axes, liquid_removal_axes];

xlim(all_axes, 'tight')
set(all_axes,'YLimitMethod','padded')

xlabel(all_layouts, "Days")
ylabel(all_layouts, "Removal efficiency [%]");
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