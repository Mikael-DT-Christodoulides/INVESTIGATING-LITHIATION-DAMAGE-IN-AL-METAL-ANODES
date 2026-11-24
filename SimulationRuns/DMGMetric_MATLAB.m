clear;
clc;

startup_mtex;

% startup commands

filename = '4.crc';

folder = 'C:\Users\mikae\OneDrive - University of Cape Town\4th Year\HonorsProject\0. FinalFiles\Data&Results\EBSD_Results\Maps';

fname = fullfile(folder, filename);

ebsd = EBSD.load(fname);


%% Gridify Data
% converts ebsd into an n x m grid (without grains), and store orientations
% and size
ebsdGrid = gridify(ebsd);    
[nRows, nCols] = size(ebsdGrid);
oriGrid = ebsdGrid.orientations;

%% Define Grains
% construct grains and store one grain id per measurement, in a separate
% variable, then make this into a grid as well
[grains,ebsdGrains] = calcGrains(ebsd,'minPixel',5);
grains = smoothBoundary(grains,5);
ebsdGrainGrid = gridify(ebsdGrains);      

%% Store grain IDs
% gId = NaN for pixels of ebsdGrid not part of any grain, otherwise the ID
% of the grain stored
gId = nan(nRows, nCols);
[isGrain, locGrain] = ismember(ebsdGrid.id, ebsdGrains.id);
gId(isGrain) = ebsdGrains.grainId(locGrain(isGrain));

% boolean grid, true wherever pixels are indexed
isIndexedGrid = reshape(ebsdGrid.isIndexed, nRows, nCols);   


 
%% Calculate GD (GOS)
% calculates overall GD value using grains from above
GD = (sum(grains.GOS .* grains.grainSize) / sum(grains.grainSize) / degree) + 1;


%% Calculate scaling factor λ and χ for each pixel
lambda = nan(nRows, nCols);
chi = nan(nRows, nCols);

for i = 1:nRows
    for j = 1:nCols
        % iteration for pixel (i, j), call it p

        pGrain = gId(i,j); 

        %* somewhat unsure of this section
        % compute N(p), cropped at map edges
        iRange = max(1,i-1):min(nRows,i+1); % x coords.
        jRange = max(1,j-1):min(nCols,j+1); % y coords.
        [iiGrid, jjGrid] = ndgrid(iRange, jRange); % row/col indices array

        neighIds = gId(iRange, jRange); % <- N(p) ∪ {p}
        neighOris = oriGrid(iRange, jRange);   % neighbour orientations
        neighIdx  = isIndexedGrid(iRange, jRange); % list of valid neighbours

        % flatten into columns, obtain p, and remove its ID from N(p), then
        % populate all other arrays
        isCentre = (iiGrid(:) == i) & (jjGrid(:) == j);
        neighIds = neighIds(:);
        neighIds(isCentre) = [];
        neighOris = neighOris(:);              
        neighOris(isCentre) = [];        
        neighIdx  = neighIdx(:);  
        neighIdx(isCentre)  = [];
    
        neighIds  = neighIds(neighIdx);        % valid ids (removed unindexed)
        neighOris = neighOris(neighIdx);       % valid orientations (removed unindexed)
        nNeighbours = numel(neighIds);         % this is now true |N(p)|
        %* 

        oriP = oriGrid(i, j); % φ(p)

        % compute χ(p)
        misSumAll = sum(angle(oriP, neighOris)) / degree;  
        chi(i,j)  = misSumAll / nNeighbours;

        % compute λ(p)

        % pixels with no gID have λ = 0
        if isnan(pGrain) || pGrain == 0
            lambda(i,j) = 0;
            continue
        end

        % determine alien pixels
        alienMask = neighIds ~= pGrain & ~isnan(neighIds) & neighIds ~= 0;   
        alienOris = neighOris(alienMask);             
        nAliens = numel(alienOris);

        % non-boundary pixels have λ = 0 as well
        if isempty(alienMask)
            lambda(i,j) = 0;
            continue
        end
        
        % since we have alien pixels, we're going to be excluding
        % some measurements (χ-λ) so we need divide by a smaller number to
        % account for this. Since we need to divide by |N(p)|-|Ap| and not 
        % |N(p)|, we'll instead just multiply χ by |N(p)|/(|N(p)|-|Ap|)

        % boundary pixel, so compute λ > 0
        misSumAlien = sum(angle(oriP, alienOris)) / degree;
        lambda(i,j) = (-(nAliens * (misSumAll - misSumAlien)) / ((nNeighbours-nAliens)*(nNeighbours)) + (misSumAlien)/(nNeighbours));
    end
end


%% Plot Figures Side by Side
figure;
mtexFig = mtexFigure;


cs = ebsd('indexed').CS;
ipfKey = ipfColorKey(cs);
ipfKey.inversePoleFigureDirection = zvector;

ipfColors = ipfKey.orientation2color(ebsd('indexed').orientations);

% IPF MAP
nextAxis;
plot(ebsd('indexed'), ipfColors, 'micronbar', 'off')
mtexTitle(sprintf('Inverse Pole Figure (IPF) Map'))

mtexFig.drawNow

% DMG Metric with KAM based on GRAINS (1) or not (2)
nextAxis;
a = 1; % scaling coefficient
dmg = (chi-(a*lambda)) * GD;
log_dmg = log10(dmg+1);

plot(ebsdGrid, log_dmg, 'micronbar', 'off')

clim([0, 3])
myMap = interp1([0 1], [0 0.39 1; 0.865 0 0], linspace(0,1,256));
colormap(gca, myMap)

mtexColorbar
mtexTitle(sprintf('DMG Metric V1.3 Map'))


mtexFig.drawNow;

disp(GD);

%outFile = fullfile(folder, 'dmg_values.txt');
%writematrix(log_dmg, outFile, 'Delimiter', 'tab');
