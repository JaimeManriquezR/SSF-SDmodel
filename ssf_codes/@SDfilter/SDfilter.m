classdef SDfilter
    % LAST MODIFIED: 23/02/2024
    properties
        domain
        epsilon
        delta
        inflow_velocity
        temperature
        light_irradiation
        mesh
    end

    methods
        function obj = SDfilter(options)
            arguments
                options.domain (1,2) {mustBeNumeric} = [-1 1];
                options.epsilon (1,1) {mustBeNumeric} = .4;
                options.delta (1,1) {mustBeNumeric} = 1e-2;
                options.inflow_velocity (1,1) {mustBeNumeric} = 2.4;
                options.temperature (1,1) {mustBeNumeric} = 293;
                options.light_irradiation (1,1) = @(t) .5*(sin(2*pi*(t - 0.25)) + 1);
                options.mesh = [];
            end

            for field = string(fieldnames(obj).')
                obj.(field) = options.(field);
            end
        end


        function obj = add_cells(obj,cell_number,options)
            arguments
                obj (1,1) SDfilter
                cell_number (1,1) {mustBeNumeric} = 0
                options.dz = [];
                options.cell_number = 1e2;
            end

            if cell_number == 0
                cell_number = options.cell_number;
            end

            if isempty(options.dz)
                dz = (obj.domain(2) - obj.domain(1))/(2*cell_number+1);
            else
                dz = options.dz;
            end
            obj.mesh.cell_centers = (obj.domain(1)+dz/2:dz:obj.domain(2)-dz/2).';
            obj.mesh.dz = dz;
        end

        % OVERRIDING FUNCTIONS
        function is_eq = isequal(self, other)
            is_eq = true;
            for field = string(fieldnames(self).')
                self_field = self.(field);
                other_field = other.(field);

                if isequal(class(self_field),'function_handle')
                    self_field = sym(self_field);
                    other_field = sym(other_field);
                end

                if ~isequal(self_field, other_field)
                    is_eq = false;
                    break;
                end
            end
        end

    end
end