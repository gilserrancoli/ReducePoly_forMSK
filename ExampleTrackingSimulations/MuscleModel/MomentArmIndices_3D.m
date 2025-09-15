% This function returns the indices of the muscles actuating the different 
% joints (for use with the moment arms) and the indices of the muscles in 
% the vector containing all muscle names.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function [Indmusi,mai] = ...
    MomentArmIndices_3D(muscleNames,muscle_spanning_joint_INFO)

NMuscle = length(muscleNames)*2;
for i = 1:length(muscleNames)
    muscleNames{i}=strrep(muscleNames{i},'_r','');
    % Muscles are ordered as left first and right second
    Indmusi.(muscleNames{i}(1:end)).l = ...
        find(strcmp(muscleNames,muscleNames{i}));
    Indmusi.(muscleNames{i}(1:end)).r = ...
        Indmusi.(muscleNames{i}(1:end)).l + NMuscle/2;
    % We select the muscles that are actuating the joints
    % In muscle_spanning_joint_INFO: the columns are organized as follows:
    % 1) hip flex, 2) hip add, 3) hip rot, 4) knee adduction, 
    % 5) knee rotation, 6) ankle flexion, 7) knee tx, 8) knee ty, 
    % 9) knee tz, 10) ankle angle, 11) subtalar, 12) lumbar ext, 
    % 13) lumbar ben, 14) lumbar rot, 15) mtp angle
    for j = 1:size(muscle_spanning_joint_INFO,2)
        if (muscle_spanning_joint_INFO(i,j) == 1)
           mai(j).mus.l(1,i) = Indmusi.(muscleNames{i}(1:end)).l;
           mai(j).mus.r(1,i) = Indmusi.(muscleNames{i}(1:end)).r;
        end        
    end    
end

for j = 1:size(muscle_spanning_joint_INFO,2)
    mai(j).mus.l(mai(j).mus.l == 0) = [];
    mai(j).mus.r(mai(j).mus.r == 0) = [];
end
end
