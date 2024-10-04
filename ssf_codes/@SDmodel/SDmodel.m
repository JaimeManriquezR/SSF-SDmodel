classdef SDmodel
    properties
        particles
        liquids
        detachment
        tau
        beta
        cahn_hilliard
        ecological
        light
    end

    methods
        function obj = SDmodel(options)
            arguments
                options.detachment = @(v,qnom) sqrt(abs(v)/qnom);
                options.tau = 1e-5;
                options.beta = .99;
                options.cahn_hilliard = struct('kappa',[],'zeta_0',[],'zeta_1',[],'mobility',[],'gradient_potential',[]);
                options.ecological = struct('stoichiometric_matrix',[],'kinetic_parameters',[],'dark_respiration',[]);
                options.light = struct('optimal',[],'attenuation',struct('particles',[],'water',[],'sand',[]));
            end
            obj.particles = cell2table(cell(0,4),'VariableNames',{'density','dispersivity','attachment_rate','transfer_rate'});
            obj.liquids = cell2table(cell(0,3),'VariableNames',{'density','dispersivity','transfer_rate'});

            obj.particles.Properties.Description = "Table containing particle components";
            obj.liquids.Properties.Description = "Table containing liquid components";

            obj.particles.Properties.VariableUnits = {'kg/m^3','m','1/day','1/day'};
            obj.liquids.Properties.VariableUnits = {'kg/m^3','m','1/day'};

            obj.particles.Properties.DimensionNames = {'Name','Parameters'};
            obj.liquids.Properties.DimensionNames = {'Name','Parameters'};

            obj.detachment = options.detachment;
            obj.tau = options.tau;
            obj.beta = options.beta;
            obj.cahn_hilliard = options.cahn_hilliard;
            obj.ecological = options.ecological;
            obj.light = options.light;
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

                if field == "cahn_hilliard"
                    for subfield = string(fieldnames(self_field).')
                        self_subfield = self_field.(subfield);
                        other_subfield = other_field.(subfield);
                        if isequal(class(self_subfield),'function_handle')
                            self_subfield = sym(self_subfield);
                            other_subfield = sym(other_subfield);
                        end
                        
                        if ~isequal(self_subfield,other_subfield)
                            is_eq = false;
                            break;
                        end

                    end

                    if is_eq
                        continue;
                    else
                        break;
                    end
                end

                if ~isequal(self_field, other_field)
                    is_eq = false;
                    break;
                end
            end
        end


    end
end