names = {'Lesi Sideboard 2/iwata_long_y100_Nx10000_results1',...
    'Lesi Sideboard 2/iwata_long_y100_Nx10000_results2_1e6',...
    'Lesi Sideboard 2/iwata_long_y100_Nx10000_results3_1e3'};

leg = {'Cutoff Size = 1e9','Cutoff Size = 1e6','Cutoff Size = 1e3'};
types = {'_delE.txt','_trec.txt'};

nn = length(names);
delEs = cell(1,nn);
trecs = cell(1,nn);

for i = 1:nn
    delEs{i} = dlmread([names{i} types{1}]);
    trecs{i} = dlmread([names{i} types{2}]);
end

fig = figure;
fig.Position = [962 42 wpic hpic];
hold on;
for i = 1:nn
    plot(delEs{i}(1:22),trecs{i}(1:22)/365,'LineWidth',2);
    hold on;
end
legend(leg,'Location','Best','FontSize',fontsizes(4));
hold off;

xlabel('Growth/Reduction Exponent Difference','FontSize',fontsizes(3));
ylabel('Recurrence Time (years)','FontSize',fontsizes(3));
axis([0.02 0.11 0 100]);
set(gca,'FontSize',fontsizes(1));
