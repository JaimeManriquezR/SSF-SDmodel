function C_sum = sum_concentrations(obj, sub_phases)
arguments
obj (1,1) SDresults
sub_phases = [{["particles","matrix"];
            ["particles", "enclosed"];
            ["liquids", "enclosed"];
            ["particles", "flowing"];
            ["liquids", "flowing"]}, {'M'; 'Pe'; 'Le'; 'Pf'; 'Lf'}]
end

C_sum = cell(length(sub_phases(:,2)),1);
[C_sum{1:end}] = deal(zeros(size(obj.frames.light)));
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
end