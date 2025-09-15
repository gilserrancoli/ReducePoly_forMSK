% This function returns the muscle indices used in the optimization problem
% as compared to all muscles of the gait2392 model
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function musi = MuscleIndices_3D_GC(muscleNames)
   
muscleNames_all =  {'addbrev','addlong','addmagProx','addmagMid',...
        'addmagDist','addmagIsch','bflh','bfsh','edl','ehl','fdl','fhl',...
        'gaslat','gasmed','gem','glmax1','glmax2','glmax3','glmed1',...
        'glmed2','glmed3','glmin1','glmin2','glmin3','grac','iliacus',...
        'pect','perbrev','perlong','pertert','piri','psoas','quadfem',...
        'recfem','sart','semimem','semiten','soleus','tfl',...
        'tibant','tibpost','vasint','vaslat',...
        'vasmed','ercspn_r','intobl_r','extobl_r',...
        'ercspn_l','intobl_l','extobl_l'};
    
count = 1;
musi = zeros(1,length(muscleNames));
for i = 1:length(muscleNames)       
    if (find(strcmp(muscleNames_all,muscleNames{i})) ~= 0)        
        musi(count) = find(strcmp(muscleNames_all,muscleNames{i}));
        count = count + 1;
    end
end

end
