function [newDefaultColors] = color_order(choice, number_of_data_sets)
% Returns a color order presets for hold all. Use these presets in the
% set function with the 'ColorOrder' options to change the default color
% order used by hold all. The list of orders is integer indexed, the key is
% as follows:
% color_order(choice, number_of_data_sets), choice variable key:
% Default => pink
%       1 => jet
%       2 => rand
%       3 => hsv
%       4 => hot
%       5 => cool
%       6 => spring
%       7 => summer
%       8 => autumn
%       9 => winter
%      10 => lines
%      11 => gray
%      12 => bone
%      13 => copper
switch choice
	case 1
		newDefaultColors = jet(number_of_data_sets);
	case 2
		newDefaultColors = rand(number_of_data_sets, 3);
	case 3
		newDefaultColors = hsv(number_of_data_sets);
	case 4
		newDefaultColors = hot(number_of_data_sets);
	case 5
		newDefaultColors = cool(number_of_data_sets);
	case 6
		newDefaultColors = spring(number_of_data_sets);
	case 7
		newDefaultColors = summer(number_of_data_sets);
	case 8
		newDefaultColors = autumn(number_of_data_sets);
	case 9
		newDefaultColors = winter(number_of_data_sets);
	case 10
		newDefaultColors = lines(number_of_data_sets);
	case 11
		newDefaultColors = gray(number_of_data_sets);
	case 12
		newDefaultColors = bone(number_of_data_sets);
	case 13
		newDefaultColors = copper(number_of_data_sets);
	otherwise
		newDefaultColors = pink(number_of_data_sets);
end
end

