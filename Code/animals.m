
function [timepoints, mask] = animals(animal_number)

if strcmp(animal_number, 'SS77')
base = '../MouseData/SS77/FunctionalActivation/Baseline';
day2 = '../MouseData/SS77/FunctionalActivation/Day2';
week1 = '../MouseData/SS77/FunctionalActivation/Week1';
week2 = '../MouseData/SS77/FunctionalActivation/Week2';
week4 = '../MouseData/SS77/FunctionalActivation/Week4';
mask = '../MouseData/SS77/brainmaskSFDI.mat';

elseif strcmp(animal_number, 'SS78')
base = '../MouseData/SS78/FunctionalActivation/Baseline';
day2 = '../MouseData/SS78/FunctionalActivation/Day2';
week1 = '../MouseData/SS78/FunctionalActivation/Week1';
week2 = '../MouseData/SS78/FunctionalActivation/Week2';
week4 = '../MouseData/SS78/FunctionalActivation/Week4';
mask = '../MouseData/SS78/brainmaskSFDI.mat';

elseif strcmp(animal_number, 'SS75')
base = '../MouseData/SS75/FunctionalActivation/Baseline';
day2 = '../MouseData/SS75/FunctionalActivation/Day2';
week1 = '../MouseData/SS75/FunctionalActivation/Week1';
week2 = '../MouseData/SS75/FunctionalActivation/Week2';
week4 = '../MouseData/SS75/FunctionalActivation/Week4';
mask = '../MouseData/SS75/brainmaskSFDI.mat';

elseif strcmp(animal_number, 'SS79')
base = '../MouseData/SS79/FunctionalActivation/Baseline';
day2 = '../MouseData/SS79/FunctionalActivation/Day2';
week1 = '../MouseData/SS79/FunctionalActivation/Week1';
week2 = '../MouseData/SS79/FunctionalActivation/Week2';
week4 = '../MouseData/SS79/FunctionalActivation/Week4';
mask = '../MouseData/SS79/brainmaskSFDI.mat';

elseif strcmp(animal_number, 'SS76')
base = '../MouseData/SS76/FunctionalActivation/Baseline';
day2 = '../MouseData/SS76/FunctionalActivation/Day2';
week1 = '../MouseData/SS76/FunctionalActivation/Week1';
week2 = '../MouseData/SS76/FunctionalActivation/Week2';
week4 = '../MouseData/SS76/FunctionalActivation/Week4';
mask = '../MouseData/SS76/brainmaskSFDI.mat';

elseif strcmp(animal_number, 'SS80')
base = '../MouseData/SS80/FunctionalActivation/Baseline';
day2 = '../MouseData/SS80/FunctionalActivation/Day2';
week1 = '../MouseData/SS80/FunctionalActivation/Week1';
week2 = '../MouseData/SS80/FunctionalActivation/Week2';
week4 = '../MouseData/SS80/FunctionalActivation/Week4';
mask = '../MouseData/SS80/brainmaskSFDI.mat';

elseif strcmp(animal_number, 'SS84')
base = '../MouseData/SS84/FunctionalActivation/Baseline';
day2 = '../MouseData/SS84/FunctionalActivation/Day2';
week1 = '../MouseData/SS84/FunctionalActivation/Week1';
week2 = '../MouseData/SS84/FunctionalActivation/Week2';
week4 = '../MouseData/SS84/FunctionalActivation/Week4';
mask = '../MouseData/SS84/brainmaskSFDI.mat';

elseif strcmp(animal_number, 'SS85')
base = '../MouseData/SS85/FunctionalActivation/Baseline';
day2 = '../MouseData/SS85/FunctionalActivation/Day2';
week1 = '../MouseData/SS85/FunctionalActivation/Week1';
week2 = '../MouseData/SS85/FunctionalActivation/Week2';
week4 = '../MouseData/SS85/FunctionalActivation/Week4';
mask = '../MouseData/SS85/brainmaskSFDI.mat';

elseif strcmp(animal_number, 'SS81')
base = '../MouseData/SS81/FunctionalActivation/Baseline';
day2 = '../MouseData/SS81/FunctionalActivation/Day2';
week1 = '../MouseData/SS81/FunctionalActivation/Week1';
week2 = '../MouseData/SS81/FunctionalActivation/Week2';
week4 = '../MouseData/SS81/FunctionalActivation/Week4';
mask = '../MouseData/SS81/brainmaskSFDI.mat';

elseif strcmp(animal_number, 'SS82')
base = '../MouseData/SS82/FunctionalActivation/Baseline';
day2 = '../MouseData/SS82/FunctionalActivation/Day2';
week1 = '../MouseData/SS82/FunctionalActivation/Week1';
week2 = '../MouseData/SS82/FunctionalActivation/Week2';
week4 = '../MouseData/SS82/FunctionalActivation/Week4';
mask = '../MouseData/SS82/brainmaskSFDI.mat';

elseif strcmp(animal_number, 'SS83')
base = '../MouseData/SS83/FunctionalActivation/Baseline';
day2 = '../MouseData/SS83/FunctionalActivation/Day2';
week1 = '../MouseData/SS83/FunctionalActivation/Week1';
week2 = '../MouseData/SS83/FunctionalActivation/Week2';
week4 = '../MouseData/SS83/FunctionalActivation/Week4';
mask = '../MouseData/SS83/brainmaskSFDI.mat';

elseif strcmp(animal_number, 'SS93')
base = '../MouseData/SS93/FunctionalActivation/Baseline';
day2 = '../MouseData/SS93/FunctionalActivation/Day2';
week1 = '../MouseData/SS93/FunctionalActivation/Week1';
week2 = '../MouseData/SS93/FunctionalActivation/Week2';
week4 = '../MouseData/SS93/FunctionalActivation/Week4';
mask = '../MouseData/SS93/brainmaskSFDI.mat';

end

timepoints = {base day2 week1 week2 week4};
end
