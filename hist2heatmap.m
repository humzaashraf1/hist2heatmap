function [Condition] = hist2heatmap(data,num,bincounts,pct)
% example data: load('data.mat')
% Requires the "surface" function in MATLAB

% Randomly samples n distributions from a large dataset and plots
% histograms for each distribution as lineplots where the colormap
% represents count intensity

% Inputs: 
    % data: 1xN cells array (row)
    % num: number of distributions to randomly sample
    % bincounts: histogram bincounts for heatmap approx
    % pct: percentile of distribution
% Outputs:
    % Condition: 3xN cell array

% Accepts a 1xN cell array where each cell represents a distribution of
% values and outputs a 3xN cell array where each row represents: 1) raw
% data values; 2) indices; and 3) coordinates for POI (nth percentile of
% each distribution)

Condition = data;
for L = 1:size(Condition,2)
    Condition{1,L} = log(Condition{1,L});
end

Condition = datasample(Condition,num,'Replace',false);
z = 1;
for L = 1:size(Condition,2)
     arsz = z.*(ones(1,size(Condition{1,L},1)))';
     sig = prctile(Condition{1,L},pct);
     Condition{2,L} = arsz;
     Condition{3,L} = {sig,Condition{2,L}(1)};
     z = z+1;
end

Condition_s = {};
cells = [];
for L = 1:num
cells(1,L) = Condition{3,L}{1,1};
end

[B,I] = sort(cells); % Sort distributions on percentile value
Condition_s(1,:) = Condition(1,I);
Condition_s(3,:) = Condition(3,I);
z = 1;
for L = 1:size(Condition_s,2)
     arsz = z.*(ones(1,size(Condition_s{1,L},1)))';
     vsig = prctile(Condition_s{1,L},pct);
     Condition_s{2,L} = arsz;
     Condition_s{3,L} = {vsig,Condition_s{2,L}(1)};
     z = z+1;
end

Condition = Condition_s;
figure()
for L = 1:size(Condition,2)
    [N,edges] = histcounts(Condition{1,L},bincounts);
    x = ones(1,bincounts)*L; z = zeros(1,bincounts);
    y = linspace(min(Condition{1,L}),max(Condition{1,L}),bincounts);
    N = rescale(N);
    surface([x;x],[y;y],[z;z],[N;N],... % Matrix interpolation for colorfill
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    colormap jet
    hold on
    scatter(Condition{3,L}{1,2},Condition{3,L}{1,1},20,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
end

vals = [];
for L = 1:size(Condition,2)
    vals(L,1) = Condition{3,L}{1};
end

%cutoff = xx;
hold on
%line([0,num],[cutoff,cutoff])
set(gca,'TickLength',[0 0])
set(gca,'FontSize',11)
ylabel('Intensity')
xlabel('Single Cells')
set(gcf,'renderer','painters') %save SVG with painters vs OpenGL render
end