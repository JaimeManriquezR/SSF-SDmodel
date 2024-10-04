classdef SDresults
    properties
        filter
        model
        frames
        time_start
        time_final
        flag
        simulation_data
    end

    methods
        function obj = SDresults(filter,model)
            arguments
                filter SDfilter = SDfilter
                model SDmodel = SDmodel
            end

            obj.filter = filter;
            obj.model = model;
            obj.frames = struct( ...
                'time', [],...
                'concentrations',[], ...
                'reactions',[], ...
                'light',[]);
            obj.time_start = 0;
            obj.time_final = 0;
            obj.flag = "UNINITIATED";
            obj.simulation_data = [];
        end

        function initial_conditions = end_results(obj)
            if isempty(obj.frames.concentrations)
                initial_conditions = [];
                return
            end

            kP = length(obj.model.particles.Name);
            kL = length(obj.model.liquids.Name);
            rhoW = obj.model.liquids.density(end);

            for j = 1:kP
                p = string(obj.model.particles.Name(j));

                initial_conditions.biofilm(:,j) = obj.frames.concentrations.matrix{p}(:,end);
                initial_conditions.biofilm(:,kP + j) = obj.frames.concentrations.enclosed{p}(:,end);
                initial_conditions.flowing(:,j) = obj.frames.concentrations.flowing{p}(:,end);
            end

            for j = 1:kL
                l = string(obj.model.liquids.Name(j));

                initial_conditions.biofilm(:,kP + kP + j) = obj.frames.concentrations.enclosed{l}(:,end);
                initial_conditions.flowing(:,kP + j) = obj.frames.concentrations.flowing{l}(:,end);
            end
            initial_conditions.enclosed_water = obj.frames.concentrations.enclosed{"WATER"}(:,end)/rhoW;
            initial_conditions.time_start = obj.time_final;

            initial_conditions.biofilm_velocity = obj.simulation_data(end).finalvelocity_biofilm;
            initial_conditions.flowing_velocity = obj.simulation_data(end).finalvelocity_flowing;
        end

        function show(obj,options)
            arguments
                obj SDresults
                options.biofilm_end = true;
                options.concentrations = false;
                options.sub_phases = false;
                options.outflow = false;
                options.light = false;
                options.biofilm_velocity = false;
                options.time_step = false;
                options.reactions = false;
                options.filtration = false;
            end

            kP = length(obj.model.particles.Name);
            kL = length(obj.model.liquids.Name);

            z = obj.filter.mesh.cell_centers;
            t = obj.frames.time;

            [T,Z] = meshgrid(t,z);
            
            subphase_type = {["particles","matrix"]; ["particles", "enclosed"]; ["liquids", "enclosed"];
                            ["particles", "flowing"];["liquids", "flowing"]};
            subphase_index = {'M'; 'Pe'; 'Le'; 'Pf'; 'Lf'};
            sub_phases = [subphase_type, subphase_index];
            
            C_sum = sum_concentrations(obj, sub_phases);

            if options.biofilm_end
                obj.plot_biofilm(concentrations = C_sum);
            end

            if options.concentrations
                obj.plot_subphases_sums(sub_phases_index = sub_phases(:,2));
            end

            if options.sub_phases
                obj.plot_subphases(sub_phases = sub_phases);
            end

            if options.outflow
                particles_index = 1:kP;
                liquids_index = 1:kL;
                obj.plot_outflow(particles_index, liquids_index);
            end

            if options.filtration
                particles_index = 1:kP;
                liquids_index = 1:kL;
                obj.plot_filtration(particles_index, liquids_index);
            end

            if options.light
                figure
                light_layout = tiledlayout(1,2);
                title(light_layout,"Light irradiation in the filter");

                nexttile;
                plot(t,obj.frames.light(1,:),'-');
                title("Incoming daylight irradiation [-]")
                xlabel("Time [d]");
                ylim([0 1])
                set(gca(),'XLimitMethod','tight')

                nexttile;
                mesh(T,Z,obj.frames.light);
                title("Effective irradiation [-]");
                xlabel("Time [d]"); ylabel("Depth [m]");
                set(gca,'YDir','reverse')
                set(gca(),'XLimitMethod','tight')
                view(2);
            end

            if options.biofilm_velocity

                tb = [];
                vb_max = [];
                fb_max = [];
                for r = 1:length(obj.simulation_data)
                    tb = [tb; obj.simulation_data(r).time(:)];
                    vb_max = [vb_max; obj.simulation_data(r).biofilm_velocity(:,1)];
                    fb_max = [fb_max; obj.simulation_data(r).biofilm_velocity(:,2)];
                end

                figure
                vel_layout = tiledlayout(2,1);

                nexttile;
                plot(tb, vb_max, 'r-', 'LineWidth', 1);
                xlabel("Time [d]");
                ylabel("Maximum biofilm velocity [m/d]");
                set(gca(),'XLimitMethod','tight')

                nexttile;
                plot(tb, fb_max, 'b-', 'LineWidth', 1);
                xlabel("Time [d]");
                ylabel("Maximum biofilm mass flux [kg/m^2d]");
                set(gca(),'XLimitMethod','tight')

                title(vel_layout, "Maximum biofilm velocity and mass flux")

            end

            if options.time_step
                iter_fig = figure;
                %dt_layout = tiledlayout(2,1);

                tb_total = [];
                plotvars_dt = {};
                plotvars_CFL = {};
                plotvars_CFLbnd = {};

                worst_CFL = [obj.simulation_data.CFL];
                worst_CFL = min([worst_CFL.worst_case]);

                it_shift = 0;
                for r = 1:length(obj.simulation_data)
                    tb = obj.simulation_data(r).time(:);
                    dt = diff(tb);
                    tb_total = [tb_total; tb(:)];

                    plotvars_dt = [plotvars_dt, {it_shift + (1:length(dt)), dt}];
                    plotvars_CFL = [plotvars_CFL, {it_shift + (1:length(dt)), obj.simulation_data(r).CFL.value}];
                    plotvars_CFLbnd = [plotvars_CFLbnd, {it_shift + [1 length(dt)], worst_CFL*[1 1]}];

                    it_shift = it_shift + length(dt);
                end

                iter_axes = axes(iter_fig);

                line_dt = semilogy(plotvars_dt{:}); hold on
                set(line_dt, ...
                    'LineStyle', '-', ...
                    'Color', 'k', ...
                    'LineWidth', 1, ...
                    'DisplayName','');

%                 line_CFL = plot(plotvars_CFL{:});
%                 set(line_CFL, ...
%                     'LineStyle', '-', ...
%                     'Color', 'b', ...
%                     'LineWidth', 1);
                line_CFLbnd = plot(plotvars_CFLbnd{:}); hold off
                set(line_CFLbnd, ...
                    'LineStyle', '--', ...
                    'Color', 'r', ...
                    'LineWidth', 1.2);

                iter_yL = ylabel(iter_axes,"$\Delta t$");
                iter_xL = xlabel(iter_axes, "N. steps");
%                 iter_L = legend([line_CFLbnd(1), line_CFL(1), line_dt(1)], ...
%                     "$\Delta t^{\rm CFL}$", ...
%                     "$\Delta t^{\rm CFL}_n$", ...
%                     "$\Delta t_{n+1}$", ...
%                     'location','east');

                iter_L = legend([line_CFLbnd(1), line_dt(1)], ...
                    "$\Delta t^{\rm CFL}$", ...
                    "$\Delta t_{n+1}$", ...
                    'location','east');


                iter_lines = iter_axes.Children;
                
                %time_axes = nexttile;
                time_fig = figure;
                time_axes = axes;
                time_lines = copyobj(iter_lines,time_axes);

                time_lines_wlegend = [];
                for l = 1:length(iter_lines)
                    time_lines(l).XData = tb_total(iter_lines(l).XData);
                    if ~isempty(iter_lines(l).DisplayName)
                        time_lines_wlegend = [time_lines_wlegend; time_lines(l)];
                    end
                end

                time_yL = ylabel(time_axes,"$\Delta t$");
                time_xL = xlabel(time_axes, "$t$");
                time_L = legend(time_lines_wlegend, ...
                    'location','east');


                time_axes.YScale = iter_axes.YScale;
                set([iter_yL, iter_xL, iter_L, time_yL, time_xL, time_L],'interpreter','latex','fontsize',14)
                set([iter_axes, time_axes],'XLimitMethod','tight')
                %title(iter_fig, "Time-stepping")
            end

            if options.reactions
                figure
                attachment_layout = tiledlayout(3,2);

                det_M = cell2mat(cellfun(@(z) sum(z,3),obj.frames.reactions{"M","attachment"},'UniformOutput',false));
                att_e = cell2mat(cellfun(@(z) sum(z,3),obj.frames.reactions{"Pe","attachment"},'UniformOutput',false));
                att_f = cell2mat(cellfun(@(z) sum(z,3),obj.frames.reactions{"Pf","attachment"},'UniformOutput',false));

                ax_e = nexttile;
                ax_f = nexttile(3);
                ax_M = nexttile(5);
                ax_net = nexttile(2,[3, 1]);

                mesh(ax_e, T, Z, att_e);
                mesh(ax_f, T, Z, att_f);
                mesh(ax_M, T, Z, det_M);
                mesh(ax_net, T,Z,att_e + att_f - det_M);

                title(ax_e, "Attachment from enclosed suspension")
                title(ax_f, "Attachment from flowing suspension")
                title(ax_M, "Detachment from biofilm matrix network")
                title(ax_net, "Net biofilm attachment/detachment")

                ax = [ax_e, ax_f, ax_M, ax_net];

                set(ax, 'XLimitMethod','tight','YDir','Reverse')
                xlabel(ax,'Time [d]')%,'interpreter','latex');
                ylabel(ax,'Depth [m]')%,'interpreter','latex');
                view(ax, 2);

                for x = ax
                    colorbar(x);
                end

                title(attachment_layout, "Attachment/detachment reactions")

                %=======================================================================%

                figure
                transfer_layout = tiledlayout(1,2);

                tr_P = cell2mat(cellfun(@(z) sum(z,3),obj.frames.reactions{"Pe","transfer"},'UniformOutput',false));
                tr_L = cell2mat(cellfun(@(z) sum(z,3),obj.frames.reactions{"Le","transfer"},'UniformOutput',false));

                ax_TP = nexttile;
                ax_TL = nexttile;

                mesh(ax_TP,T,Z,tr_P);
                mesh(ax_TL,T,Z,tr_L);

                title(ax_TP, "Transfer of particlesfrom flowing to enclosed")
                title(ax_TL, "Transfer of liquids from flowing to enclosed")

                ax = [ax_TP, ax_TL];
                set(ax, 'XLimitMethod','tight','YDir','Reverse')
                xlabel(ax,'Time [d]')%,'interpreter','latex');
                ylabel(ax,'Depth [m]')%,'interpreter','latex');
                view(ax, 2);

                for x = ax
                    colorbar(x);
                end

                title(transfer_layout, "Transfer reactions")

                %=======================================================================%

                for i = 1:kP
                    figure
                    ecological_layout(i) = tiledlayout(3,2);

                    rx = cellfun(@(z) z(:,:,i), obj.frames.reactions{["M"; "Pe"; "Pf"],"ecological"}, 'UniformOutput', false);

                    ax_b = nexttile;
                    ax_e = nexttile(3);
                    ax_f = nexttile(5);
                    ax_T = nexttile(2,[3, 1]);

                    mesh(ax_b, T, Z, rx{1});
                    mesh(ax_e, T, Z, rx{2});
                    mesh(ax_f, T, Z, rx{3});
                    mesh(ax_T, T, Z, sum(cat(3,rx{:}),3));

                    title(ax_b, "Biofilm matrix network");
                    title(ax_e, "Enclosed suspension");
                    title(ax_f, "Flowing suspension")
                    title(ax_T, "All subphases");
                    title(ecological_layout(i),"Ecological reactions");
                    subtitle(ecological_layout(i),obj.model.particles.Name{i});

                    ax = [ax_b, ax_e, ax_f, ax_T];

                    set(ax, 'XLimitMethod','tight','YDir','Reverse')
                    xlabel(ax,'Time [d]')%,'interpreter','latex');
                    ylabel(ax,'Depth [m]')%,'interpreter','latex');
                    view(ax,2);

                    for x = ax
                        colorbar(x);
                    end

                end

            end


        end

        function [phi,z,t] = volume_fractions(obj)
            rhoP = obj.model.particles.density(1);
            rhoL = obj.model.liquids.density(1);

            z = obj.filter.mesh.cell_centers;
            t = obj.frames.time;
            
            subphase_type = {["particles","matrix"]; ["particles", "enclosed"]; ["liquids", "enclosed"];
                            ["particles", "flowing"];["liquids", "flowing"]};
            subphase_index = {'M'; 'Pe'; 'Le'; 'Pf'; 'Lf'};
            sub_phases = [subphase_type, subphase_index];

            C_sum = cell(length(sub_phases(:,2)),1);
            [C_sum{1:end}] = deal(zeros(length(z),length(t)));
            C_sum = table(C_sum{:}, ...
                'VariableNames',sub_phases(:,2));
            for s = 1:size(sub_phases,1)
                component_type = sub_phases{s,1}(1);
                volume_type = sub_phases{s,1}(2);
                phase_name = sub_phases{s,2};

                for i = 1:length(obj.model.(component_type).Name)
                    C_sum.(phase_name) = C_sum.(phase_name) + ...
                        obj.frames.concentrations{obj.model.(component_type).Name{i},volume_type}{:};
                end
            end

            phiM = C_sum.("M")/rhoP;
            phie = C_sum.("Pe")/rhoP + C_sum.("Le")/rhoL;
            phib = phiM + phie;

            phi = struct('biofilm',phib, ...
                'enclosed_suspension',phie, ...
                'biofilm_matrix',phiM);
        end


        function obj = concatenate(obj,new_obj)
            if isequal(obj.flag,"UNINITIATED")
                obj = new_obj;
                return
            end

            if isequal(new_obj.flag,"UNINITIATED")
                return
            end
           

            if isequal(obj.filter,new_obj.filter) && isequal(obj.model, new_obj.model)
                if isequal(obj.time_final, new_obj.time_start)
                    obj.frames.time = [obj.frames.time; new_obj.frames.time(2:end)];
                    for field = ["concentrations", "reactions"]
                        for i = 1:size(obj.frames.(field),1)
                            for j = 1:size(obj.frames.(field),2)
                                if isequal(obj.frames.(field){i,j}, {0})
                                    continue;
                                end
                                obj.frames.(field){i,j} = ...
                                    {[obj.frames.(field){i,j}{:}, new_obj.frames.(field){i,j}{:}(:,2:end,:)]};
                            end
                        end
                    end
                    obj.frames.light = [obj.frames.light, new_obj.frames.light(:,2:end)];

                    obj.time_final = new_obj.time_final;
                    obj.flag = new_obj.flag;

                    obj.simulation_data = [obj.simulation_data, new_obj.simulation_data];

                else
                    error("Initial time of concatenated object does not coincide with final time of the original.")
                end
            else
                error("Filter or model not compatible between results objects.")
            end
        end

        function [new_filter, inoculated_conditions, scraped_m] = scrape(obj,z_s,geometry_preserve)
            arguments
                obj SDresults
                z_s (1,1) {mustBeNumeric} = 5e-2 % default: 5 cm
                geometry_preserve = false
            end
            end_conditions = obj.end_results;
            old_filter = obj.filter;

            dz = old_filter.mesh.dz;
            scraped_cells = ceil(z_s/dz);
            scraped_m = scraped_cells*dz;

            new_domain = old_filter.domain - scraped_m;
            new_mesh.cell_centers = old_filter.mesh.cell_centers - scraped_m;
            new_mesh.dz = old_filter.mesh.dz;

            new_filter = old_filter;
            if ~geometry_preserve
                new_filter.domain = new_domain;
                new_filter.mesh = new_mesh;
            end

            inoculated_conditions = end_conditions;

            % clean flowing suspension
            inoculated_conditions.time_start = 0;
            inoculated_conditions.biofilm_velocity = [];
            inoculated_conditions.flowing_velocity = [];

            % scrape above sand
            supernatant_water = new_mesh.cell_centers <= dz/2;
            inoculated_conditions.biofilm(supernatant_water,:) = 0;
            inoculated_conditions.enclosed_water(supernatant_water,:) = 0;
            inoculated_conditions.flowing(supernatant_water,:) = repmat(obj.simulation_data(end).inflow_concentrations(0),[nnz(supernatant_water) 1]);
        end


    end
end