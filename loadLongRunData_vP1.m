%% loadLongRunData_v1_0
%  version 1.0
%  Author: Adeyinka Lesi
%  Date: 8/7/18
%  for long run with treatment only

fileList = {...
    'iwata_long_y100_Nx10000_10172018-1006_1',...
    'iwata_long_y100_Nx10000_10172018-1006_2'};
    
bp = [0.64802,0.80825,0,0,4.602e-08,4.602e-11,0.91173,Inf,0,0];
rec_thres = 1e12;

delEs = zeros(1,length(fileList)); 
tfin = 80;
lims = [0 tfin 0 4e12];
    
afig = figure;
afig.Position = [962 42 wpic hpic];

for i = 1:length(fileList)
    load(fileList{i});
    res1 = mg_res;
    key1 = res1.key;
    t3 = res1.t3;
    t2 = res1.t2;
    delEs(i) = bp(2)-key1.DEATH_PARAMETER2;
    % total number of cells
    itfin = round((tfin*365-t3(1))/(t3(2)-t3(1)))+1;
    plot(res1.t3(1:itfin)/365,res1.num_cells(1:itfin),'LineWidth',2);
    hold on;
end
plot(res1.t3(1:itfin)/365,ones(1,itfin)*rec_thres,'k--');

axis(lims);
xlabel('Time (years)','FontSize',fontsizes(2));
ylabel('Total Number of Tumor Cells','FontSize',fontsizes(3))
text(4,2e12,'Surgery after 3 years','FontSize',fontsizes(4));

leg = cell(1,length(fileList)+1);
for j = 1:length(fileList)
    leg{j} = sprintf('Exponent Diff = %3.2f',delEs(j));
end
leg{end} = 'Recurrence Threshold';
legend(leg,'Location',[0.2 0.75 0.5 0.15],'FontSize',fontsizes(4));
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
legend boxoff;



    