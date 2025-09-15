% This function returns the GRFs given the path to the GRF file
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function GRF = getGRF(pathGRF)

load(pathGRF,'GRFall')
GRF.time = GRFall.data(:,contains(GRFall.colheaders,{'time'}));
if contains(pathGRF,'_tm_')
    identifier ={'grforce1_f','grforce2_f'};
    identifierp={'grforce1_p','grforce2_p'};
    identifierm={'grforce1_t','grforce2_t'};
else
    identifier = {'grforce2_f','grforce1_f','grforce3_f'};     % Right leg: grforce2; Left leg: grforce1; Right leg: grforce3
    identifierp = {'grforce2_p','grforce1_p','grforce3_p'};    % Right leg: grforce2; Left leg: grforce1; Right leg: grforce3
    identifierm = {'grforce2_t','grforce1_t','grforce3_t'};    % Right leg: grforce2; Left leg: grforce1; Right leg: grforce3
end
GRF.val.all(:,1) = GRF.time;
GRF.pos.all(:,1) = GRF.time;
GRF.Mcop.all(:,1) = GRF.time;
GRF.MorGF.all(:,1) = GRF.time;

GRF.val.all3(:,1) = GRF.time;
GRF.pos.all3(:,1) = GRF.time;
GRF.Mcop.all3(:,1) = GRF.time;
GRF.MorGF.all3(:,1) = GRF.time;
axis = {'x','y','z'};
if contains(pathGRF,'_tm_')
    leg={'l','r1'};
else
    leg = {'r1','l','r2'};
end
count = 1;
for j = 1:length(leg)
    for i = 1:length(axis)    
        count = count + 1;
        GRF.val.(leg{j})(:,i) = GRFall.data(:,contains(GRFall.colheaders,[identifier{j},axis{i}]));
        GRF.val.all3(:,count) = GRF.val.(leg{j})(:,i);
        GRF.pos.(leg{j})(:,i) = GRFall.data(:,contains(GRFall.colheaders,[identifierp{j},axis{i}]));
        GRF.pos.all3(:,count) = GRF.pos.(leg{j})(:,i);
        GRF.Mcop.(leg{j})(:,i) = GRFall.data(:,contains(GRFall.colheaders,[identifierm{j},axis{i}]));
        GRF.Mcop.all3(:,count) = GRF.Mcop.(leg{j})(:,i);
    end
end
intervr1=find(GRF.val.all3(:,3)>10);
GRF.val.all(intervr1,2:4)=GRF.val.all3(intervr1,2:4);
GRF.val.all(:,5:7)=GRF.val.all3(:,5:7);
if contains(pathGRF,'_tm_')
    GRF.Mcop.Y3=GRF.Mcop.all3(:,[3,6]);
else
    intervr2=find(GRF.val.all3(:,9)>10);
    GRF.val.all(intervr2,2:4)=GRF.val.all3(intervr2,8:10);
    GRF.Mcop.Y3 = GRF.Mcop.all3(:,[3,6,9]); 
end
GRF.Mcop.Y(:,2)=GRF.Mcop.Y3(:,2);
GRF.Mcop.Y(intervr1,1)=GRF.Mcop.Y3(intervr1,1);
if contains(pathGRF,'_tm_')
else
    GRF.Mcop.Y(intervr2,1)=GRF.Mcop.Y3(intervr2,3);
end

% Calculate moments with respect to ground frame origin
for j = 1:length(leg)
    GRF.MorGF.(leg{j})(:,1) =  GRF.pos.(leg{j})(:,2).*GRF.val.(leg{j})(:,3) - GRF.pos.(leg{j})(:,3).*GRF.val.(leg{j})(:,2);
    GRF.MorGF.(leg{j})(:,2) =  GRF.pos.(leg{j})(:,3).*GRF.val.(leg{j})(:,1) - GRF.pos.(leg{j})(:,1).*GRF.val.(leg{j})(:,3) + GRF.Mcop.(leg{j})(:,2);
    GRF.MorGF.(leg{j})(:,3) =  GRF.pos.(leg{j})(:,1).*GRF.val.(leg{j})(:,2) - GRF.pos.(leg{j})(:,2).*GRF.val.(leg{j})(:,1);
end
GRF.MorGF.all3(:,2:4) = GRF.MorGF.r1;
GRF.MorGF.all3(:,5:7) = GRF.MorGF.l;
if contains(pathGRF,'_tm_')
else
    GRF.MorGF.all3(:,8:10) = GRF.MorGF.r2;
end

GRF.MorGF.all(:,5:7)=GRF.MorGF.all3(:,5:7);
GRF.MorGF.all(intervr1,2:4)=GRF.MorGF.all3(intervr1,2:4);
if contains(pathGRF,'_tm_')
else
    GRF.MorGF.all(intervr2,2:4)=GRF.MorGF.all3(intervr2,8:10);
end
end
