clc
clear variables

%% Part 1: Preparation
% First I get a list of all the .tiff files in the folder and initiate a
% matrix that has a row for each image analysis.
% Function dir() gives me a struct and fullfile() is used to build file 
% paths (see MATLAB documentation)

folderPath = 'C:\Users\miri\OneDrive\Desktop\MATLAB\Images\CDNB15';

fileList = dir( ...
    fullfile( ...
    folderPath, '*.tiff')); 

numImages = numel(fileList); % for image iteration

% Initiate matrix to store the results 
% (column value = #parameters to be saved)
stats = zeros(numImages, 3); 

%% Part 2: Image segmentation
% Next I loop through each of the files and process them one by one.
% .name — goes to field File name in struct
for i = 1:numImages
    currentFile = imread(fullfile(folderPath, fileList(i).name)); 

    % Colors increase complexity of the image and makes image segmentation
    % a lot harder. After conversion to grayscale I adjust the contrast to
    % enhance blurry cell edges.
    grayImage = im2gray(currentFile);

    enhancedImage = adapthisteq( ...
        grayImage, ...
        'ClipLimit', 0.007); 

    % Function imbinarize() obtains a binary image; automatic threshold 
    % detection with graythresh()
    binaryImage = imbinarize( ...
        enhancedImage, ...
        graythresh(enhancedImage));

    % Edge detection and fill the holes 
    % --> obtains bw (background black, cells white)
    edgeImage = edge(binaryImage, 'canny');
    filledImage = imfill(edgeImage, 'holes');

    % Erosion and dilation with imopen() to get rid of small objects. 
    % For strel n = 0: the structuring element members comprise all pixels 
    % whose centers are no greater 
    % than r away from the origin - instead of approximation.
    cells = strel('disk', 40, 0); 
    filteredImage = imopen(filledImage, cells);

    % Use watershed algorithm to separate the cells --> IMPORTANT for the
    % component analysis because afterwarts the ridgeline between the cells
    % will be set to zero. Therefore no connection will be detected with
    % bwconncomp() function.
    bw1 = -bwdist(~filteredImage);
    wsImage = watershed(bw1);
    bw = filteredImage;
    bw(wsImage == 0) = 0;
    mask = imextendedmin(bw1,2);

    bw2 = imimposemin(bw1,mask);
    ws2_control = watershed(bw2);
    bw3 = filteredImage;
    bw3(ws2_control == 0) = 0;
 
    % Another round of erosion and dilation is necessary because otherwise
    % I have the ridgelines beween the cells that change their shape.
    bwRBC = imopen(bw3, cells);
    erodedImage = imerode(bw3, strel('disk', 20, 0));

    
%% Part 3: Feature extraction
% Here I can use the previously created mask (maskRBC) for the analysis
    separatedCC1 = bwconncomp(bwRBC);
    separatedCC2 = bwconncomp(erodedImage);
    numSeparatedCells = separatedCC1.NumObjects;

    % Array to store information about individual cells
    cellInfo = zeros(numSeparatedCells, 3);
    
    props = regionprops(separatedCC1, 'Circularity');
    circularity = struct2array(props);

    numStomatocytic = 0; 
    numBiconcave = 0;

    % Display the original grayscale image with separated cells outlined
    figure;
    imshow(currentFile, []);

    hold on

        % Loop through the separated cells extract desired information 
        % and outline them
        for j = 1:numSeparatedCells
        cellInfo(j, 1) = j; 
        cellInfo(j, 3) = circularity(1, j); 
        cellPixels1 = separatedCC1.PixelIdxList{j};
        cellPixels2 = separatedCC2.PixelIdxList{j};
    
        % Extract pixel values for the current cell from the grayscale 
        % image
        cellValues2 = grayImage(cellPixels2);
        meanValues = mean(cellValues2);
        cellInfo(j, 2) = meanValues;
    
        % Calculate the coefficient of variation (CV) the cells pixel 
        % intensities
        CV = std(double(cellValues2))/meanValues;
        cellInfo(j, 4) = CV;
  
        [rows, cols] = ind2sub(size(bwRBC), cellPixels1);
        boundary = bwtraceboundary(bwRBC, [rows(1), cols(1)], 'N');

            if cellInfo(j, 3) > 0.94 && cellInfo(j, 4) >= 0.07
                plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
                numStomatocytic = numStomatocytic + 1;
            elseif cellInfo(j, 3) > 0.97 && cellInfo(j, 4) <= 0.045
                plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2);
                numBiconcave = numBiconcave + 1;
            else 
                plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 2);
            end

        % Display the CV near the cell
        numNone = numSeparatedCells - (numStomatocytic + numBiconcave);
        stats(i, 1) = numStomatocytic;
        stats(i, 2) = numBiconcave;
        stats(i, 3) = numNone;
        text(cols(1), rows(1), sprintf('CV: %.2f', cellInfo(j, 4)), ...
            'Color', 'b', 'FontSize', 8, 'FontWeight', 'normal');
        end

    hold off 
end

%% Part 4: Save the results in an .xlsx file
% Because I want to have column names I decided that the easiest way is 
% to create a table from stats 
columnNames = {'Stomatocytic', 'Biconcave', 'k. A.'};
t = array2table(stats, 'VariableNames', columnNames);

writetable(t, 'C:\Users\miri\OneDrive\Desktop\MATLAB\Images\analysis.xlsx');

