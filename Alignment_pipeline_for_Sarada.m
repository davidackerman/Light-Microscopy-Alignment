%% Parameters that are more likely to be changed
animalId = '382366';
scope = 'pe'; %pe
stackProject = 'looger_mouse';
orderedInputDir = ['/groups/looger/loogerlab/render/data/multi_channel/animal/' animalId '/'];
baseOutputDir =  '/groups/looger/loogerlab/render/data/split_channel/';%'/groups/flyTEM/home/ackermand/looger_lab/render/data/split_channel/'; %
nfirst = 1;
nlast = Inf; %If nlast = Inf, it will run up to and including the last image
%% Parameters that are unlikely to change
alignmentDiagnosticsDirFine = ['/groups/looger/loogerlab/render/data/alignment_diagnostics/fine/animal/' animalId '/'];%['/groups/flyTEM/home/ackermand/looger_lab/render/data/alignment_diagnostics/fine/animal/' animalId '/'];%
alignmentDiagnosticsDirFiner = ['/groups/looger/loogerlab/render/data/alignment_diagnostics/finer/animal/' animalId '/'];%['/groups/flyTEM/home/ackermand/looger_lab/render/data/alignment_diagnostics/fine/animal' animalId '/'];%
scapeDirectory = '/groups/looger/loogerlab/render/scapes/';%'/groups/flyTEM/home/ackermand/looger_lab/render/scapes/';%'/groups/looger/loogerlab/render/scapes/';%
if ~strcmp(baseOutputDir(end),'/'), baseOutputDir(end+1) = '/'; end
if ~strcmp(scapeDirectory(end),'/'), scapeDirectory(end+1) = '/'; end
timing = tic;
currentDirectory = pwd;
numPixelsScaleFine = 200000;

channels = {'dapi', 'red', 'green'};
channelOrder = [3,1,2]; %eg R, G, B

stackOwner = 'flyTEM';
stackResX = 50.0;
stackResY = 50.0;
stackResZ = 100.0;
renderHostAndPort = 'tem-services.int.janelia.org:8080';
baseDataUrl = ['http://' renderHostAndPort '/render-ws/v1'];
baseStackName = [animalId '_acquire'];

rc.stack = [baseStackName '_dapi'];
rc.owner = stackOwner;
rc.project = stackProject;
rc.service_host = '10.40.3.162:8080';
rc.baseURL        = ['http://' rc.service_host '/render-ws/v1'];
rc.verbose = 0;

if isequal(nlast, inf), nlast = numel(dir([orderedInputDir '*.tif'])); end
    
rcfine = rc;
rcfine.stack = [animalId '_dapi_fine_' num2str(nfirst) '_' num2str(nlast)];

rcfiner = rcfine;
rcfiner.stack = [animalId '_dapi_finer_' num2str(nfirst) '_' num2str(nlast)];

pmfine.server           = 'http://10.40.3.162:8080/render-ws/v1';
pmfine.owner            = 'loogerl';%'loogerl'; %change to loogerl
pmfine.match_collection = ['mouse_' animalId '_dapi_fine_' num2str(nfirst) '_' num2str(nlast)];
pmfine.npoints_per_pair = 200;

pmfiner = pmfine;
pmfiner.match_collection = ['mouse_' animalId '_dapi_finer_' num2str(nfirst) '_' num2str(nlast)];

tile.server = pmfine.server;
tile.owner = rcfiner.owner;
tile.project = rcfiner.project;
tile.stack = rcfiner.stack;

visibleAlignmentDiagnosticFigures = false;

outputScale = 1.0;
%% Create stack and do alignment
%Create collection

if baseOutputDir(end)~='/', baseOutputDir(end+1) = '/'; end
if orderedInputDir(end)~='/', orderedInputDir(end+1) = '/'; end

% Create necessary directories
resp = system(['mkdir -p ' alignmentDiagnosticsDirFine]);
resp = system(['mkdir -p ' alignmentDiagnosticsDirFiner]);
splitChannelOutputDir = [baseOutputDir  scope  '/' animalId '/'];
tileSpecsOutputDir = [splitChannelOutputDir 'tileSpecs/'];
resp = system(['mkdir -p ' tileSpecsOutputDir]);

% Create a collection for each channel
parfor channelIndex = 1:3
    rccurrent = rc;
    currentStackName = getChannelStackName(baseStackName, channels{channelIndex});
    rccurrent.stack = currentStackName;
    delete_renderer_stack(rccurrent);
    createStack(baseDataUrl, stackOwner, stackProject, currentStackName, stackResX, stackResY, stackResZ);
end

% Split images and create array of tile specs
normalizeImageNames(orderedInputDir);
[status, cmdout] = system(['ls ' orderedInputDir '*.tif']);
if status~=0, error(['No images found in directory ' orderedInputDir]); end
tifImageNames = strsplit(cmdout,'\n');
tifImageNames(end) = [];

tileSpecsArray = cell(3, numel(tifImageNames));
fprintf(['\nImage splitting progress:\n' repmat('.',1,numel(tifImageNames)) '\n\n']);
imageSizes = zeros(numel(tifImageNames),3);
originalImageTypes = cell(numel(imageSizes),1);
parfor z = 1:numel(tifImageNames)
    path = tifImageNames{z};
    [tileSpecsForEachChannel, imageHeight, imageWidth, originalImageTypes{z}] = splitMultiChannelImage(path, splitChannelOutputDir, z);
    imageSizes(z,:) = [imageWidth, imageHeight, imageWidth*imageHeight];
    tileSpecsArray(:,z) = tileSpecsForEachChannel;
    fprintf('\b|\n');
end
% Put restrictions on the minimum scale
originalImageType = originalImageTypes{1};
minimumImageDimension = min(min(imageSizes(:,1:2)));
minimumScale = 256/minimumImageDimension;
averageImageSize = mean(imageSizes(:,3));
scaleFine = max(min(1,sqrt(numPixelsScaleFine/averageImageSize)),minimumScale);
scaleFiner = min(1,2*scaleFine);
rcfine.versionNotes = ['Fine alignment using matrix solver with Matlab backslash operator -- fast script. Scale ' num2str(scaleFine) '. Alignment diagnostic figures located in ' alignmentDiagnosticsDirFine];
rcfiner.versionNotes = ['Finer alignment using matrix solver with Matlab backslash operator -- fast script. Scale ' num2str(scaleFiner) '. Alignment diagnostic figures located in ' alignmentDiagnosticsDirFine];

% Write the tile specs to corresponding json files
writeTileSpecsJsonFile([tileSpecsOutputDir 'red_sections.json'], tileSpecsArray(1,:));
writeTileSpecsJsonFile([tileSpecsOutputDir 'green_sections.json'], tileSpecsArray(2,:));
writeTileSpecsJsonFile([tileSpecsOutputDir 'dapi_sections.json'], tileSpecsArray(3,:));

% Import the json specs and complete the stack
rccurrent = rc;
for channelIndex = 1:3
    currentStackName = getChannelStackName(baseStackName, channels{channelIndex});
    rccurrent.stack = currentStackName;
    jsonFile = [tileSpecsOutputDir channels{channelIndex} '_sections.json'];
    importJson(jsonFile, baseDataUrl, stackOwner, stackProject, currentStackName);
    set_renderer_stack_state_complete(rccurrent);
end

% Do fine alignments using two different scales
transfac = 10^-12;
lambda_exponent=7;
visibility = 'off';
fineAlign(nfirst, nlast, rc, rcfine, pmfine, scaleFine, transfac, lambda_exponent, alignmentDiagnosticsDirFine,'fine',visibility);
fineAlign(nfirst, nlast, rcfine, rcfiner, pmfiner, scaleFiner, transfac, lambda_exponent, alignmentDiagnosticsDirFiner,'finer',visibility);

% Create the other fine aligned channels by copying the tile specs from the
% finer aligned channel
for channelIndex = 2:3
    color = channels{channelIndex};
    currentStackName = [animalId '_' color '_finer_' num2str(nfirst) '_' num2str(nlast)];
    rccurrent = rcfiner;
    rccurrent.stack = currentStackName;
    delete_renderer_stack(rccurrent);
    createStack(baseDataUrl, stackOwner, stackProject, currentStackName, stackResX, stackResY, stackResZ);
    channelTileSpecsArray = cell(1, nlast-nfirst+1);
    count = 1;
    for z=nfirst:nlast
        [~, imgName, ~] = fileparts(tifImageNames{z});
        tile.renderer_id = [ imgName '.' num2str(z) '.0'];
        currentChannelRendererID = [ imgName '.' num2str(z) '.0'];
        dapiTileSpec = get_tile_spec_renderer(tile);
        imagePath = strrep(strrep(dapiTileSpec.tileSpecs.mipmapLevels.x0.imageUrl, 'file:',''), 'dapi', color);
        copyFromDAPI = true;
        channelTileSpecsArray{count} = generateTileSpec(currentChannelRendererID, z, imagePath, dapiTileSpec.tileSpecs.width, dapiTileSpec.tileSpecs.height, dapiTileSpec); 
        count=count+1;
    end
    finerJsonFile = [tileSpecsOutputDir color '_finer_sections.json'];
    writeTileSpecsJsonFile(finerJsonFile, channelTileSpecsArray)
    importJson(finerJsonFile, baseDataUrl, stackOwner, stackProject, currentStackName);
    set_renderer_stack_state_complete(rccurrent);
end

createOutputImages( rcfiner, channels, channelOrder, scapeDirectory, outputScale, originalImageType)

cd(currentDirectory);
totalTime = toc(timing);
%% Functions

function channelStackName = getChannelStackName(baseStackName, channelName)
% Make the channel stack name
channelStackName = [baseStackName '_' channelName];
end

function resp = createStack(baseDataUrl, stackOwner, stackProject, currentStackName, stackResX, stackResY, stackResZ)
systemCommand = sprintf('/groups/flyTEM/flyTEM/render/bin/manage_stacks.sh --action CREATE --baseDataUrl %s --owner %s --project %s --stack %s --stackResolutionX %0.1f --stackResolutionY %0.1f --stackResolutionZ %0.1f',...
    baseDataUrl, stackOwner, stackProject, currentStackName, stackResX, stackResY, stackResZ);
resp = system(systemCommand);
end

function normalizeImageNames(orderedInputDir)
% Normalize the image names
tifImageNames = dir([orderedInputDir '*.tif']);
parfor i=1:numel(tifImageNames)
    [~,unnormalizedName,ext] = fileparts(tifImageNames(i).name);
    unnormalizedNameForMoving = regexprep(unnormalizedName,' ','\\ ');
    normalizedName = zeroPadNumbersInString(regexprep(unnormalizedName,'[^a-zA-Z0-9_]','_'));
    [~,~] = system(['mv ' orderedInputDir unnormalizedNameForMoving ext ' ' orderedInputDir normalizedName ext]);
end
end

function outputString  = zeroPadNumbersInString(inputString)
% Make sure numbers are zero padded
numberIndices = regexp(inputString,'\d');
outputString = [];
for i = numel(inputString):-1:1
    outputString = [inputString(i), outputString];
    if ismember(i, numberIndices) 
        if ~ismember(i-1, numberIndices) %then might need to add zeros
            if ~ismember(i+1, numberIndices) %if the next value is not a number, need to pad with two zeros
                outputString = ['00', outputString];
            elseif ~ismember(i+2, numberIndices)
                outputString = ['0', outputString];
            end
        end
    end
end
end

function  importJson(jsonFile, baseDataUrl, stackOwner, stackProject, currentStackName)
% Import json file
    systemCommand = sprintf(['/groups/flyTEM/flyTEM/render/bin/import_json.sh %s --baseDataUrl %s --owner %s --project %s --stack %s '...
        '--validatorClass org.janelia.alignment.spec.validator.TemTileSpecValidator --validatorData minCoordinate:-500000,maxCoordinate:500000,minSize:1,maxSize:500000'],...
        jsonFile, baseDataUrl, stackOwner, stackProject, currentStackName);
    [~, cmdoutput] = system(systemCommand);
    if strfind(cmdoutput, 'ERROR'), error(['Error in importing json file ' jsonFile]); end
end

function tileSpecJson = saveSplitImage(img, commonSplitPath, channelName, imgName, z, imgWidth, imgHeight)
% Saves the split image and create4s the corresponding tile spec
splitImagePath = [commonSplitPath '__' channelName '.tif'];
imwrite(img, splitImagePath, 'Compression', 'none');

tileId = [imgName '.' num2str(z) '.0'];
tileSpecJson = generateTileSpec(tileId, z, splitImagePath, imgWidth, imgHeight);
end

function [tileSpecsForEachChannel, imgHeight, imgWidth, imageClass] = splitMultiChannelImage(multiChannelImagePath, splitChannelOutputDir, z)
% Splits the image and stores the tile spec for each channel in a cell
% array
[~,imgName,~]= fileparts(multiChannelImagePath);
commonSplitPath = [splitChannelOutputDir imgName];
img = imread(multiChannelImagePath);
[imgHeight, imgWidth] = size(img(:,:,1)); 

tileSpecsForEachChannel = cell(3,1);
tileSpecsForEachChannel{3} = saveSplitImage(img(:,:,3), commonSplitPath, 'dapi', imgName, z, imgWidth, imgHeight);
tileSpecsForEachChannel{1} = saveSplitImage(img(:,:,1), commonSplitPath, 'red', imgName, z, imgWidth, imgHeight);
tileSpecsForEachChannel{2} = saveSplitImage(img(:,:,2), commonSplitPath, 'green', imgName, z, imgWidth, imgHeight);

imageClass = class(img);
end

function tileSpecString = generateTileSpec(tileID, z, imagePath, width, height, dapiTileSpec)
if nargin < 6
    dataString = '1 0 0 1 0 0'; 
else
    dataString = dapiTileSpec.tileSpecs.transforms.specList.dataString;
end
% Creates the tile spec for given input
% 1: tileId, 2: z, 3: imagePath, 4: width, 5: height
tileSpecString = sprintf([
    '{\n'...
    '  "tileId" : "%s",\n'...
    '  "layout" : {\n'...
    '    "sectionId" : "%d.0",\n'...
    '    "temca" : "0",\n'...
    '    "camera" : "0",\n'...
    '    "imageRow" : 0,\n'...
    '    "imageCol" : 0,\n'...
    '    "stageX" : 0.0,\n'...
    '    "stageY" : 0.0,\n'...
    '    "rotation" : 0.0\n'...
    '  },\n'...
    '  "z" : %d.0,\n'], tileID, z, z);
if nargin==6
    fromDapi = sprintf([
        '  "minX" : %d.0,\n' ...s%s/%s/scale_%0.1f_align/000/
        '  "minY" : %d.0,\n' ...
        '  "maxX" : %d.0,\n' ...
        '  "maxY" : %d.0,\n' ...
        ], dapiTileSpec.tileSpecs.minX, dapiTileSpec.tileSpecs.minY, dapiTileSpec.tileSpecs.maxX, dapiTileSpec.tileSpecs.maxY);
    tileSpecString = [tileSpecString fromDapi];
end
    lastChunk = sprintf([
    '  "width" : %d.0,\n'...
    '  "height" : %d.0,\n'...
    '  "mipmapLevels" : {\n'...
    '    "0" : {\n'...
    '      "imageUrl" : "file:%s"\n'...
    '    }\n'...
    '  },\n'...
    '  "transforms" : {\n'...
    '    "type" : "list",\n'...
    '    "specList" : [ {\n'...
    '      "type" : "leaf",\n'...
    '      "className" : "mpicbg.trakem2.transform.AffineModel2D",\n'...
    '      "dataString" : "%s"\n'...
    '    } ]\n'...
    '  }\n'...
    '}'], width, height, imagePath, dataString);
tileSpecString = [tileSpecString lastChunk];
end


function writeTileSpecsJsonFile(jsonFilePath, tileSpecs)
% Writes the tile specs to a json file
fileID = fopen(jsonFilePath, 'w');
fprintf(fileID,'[\n');
for tileNumber=1:numel(tileSpecs)
    fprintf(fileID, tileSpecs{tileNumber});
    if tileNumber~=numel(tileSpecs), fprintf(fileID, ','); end
    fprintf(fileID, '\n');
end
fprintf(fileID, ']\n');
fclose(fileID);
end

% Fine alignment
function fineAlign(nfirst, nlast, rc, rcfine, pmfine, scale, transfac, lambda_exponent, alignmentDiagnosticsDir, fineOrFiner, visibility)
diary on; clc;

pm_fine_scale = scale; %0.1250;%0.25
dir_scratch = ['/scratch/' char(java.lang.System.getProperty('user.name')) '/'];
% % ------------------------------------------------------------------------------------
 
% configure solver
opts.degree = 1;    % 1 = affine, 2 = second order polynomial, maximum is 3
opts.solver = 'backslash';% %'pastix';%%'gmres';%'backslash';'pastix';

opts.transfac = transfac;  % translation parameter regidity
opts.xfac = 1;   % 2nd order parameter rigidity in x
opts.yfac = 1;   % 2nd order parameter regidity in y
opts.lambda = 10^(lambda_exponent); % 10^4.5 (best results so far) ------------------>

opts.nbrs = 1;
opts.nbrs_step = 1;
opts.xs_weight = 1.0;
opts.min_points = 5;
opts.max_points = inf;
opts.filter_point_matches = 0;

opts.outlier_lambda = 1e2;  % large numbers result in fewer tiles excluded
opts.min_tiles = 2; % minimum number of tiles that constitute a cluster to be solved. Below this, no modification happens
opts.matrix_only = 0;   % 0 = solve , 1 = only generate the matrix
opts.distribute_A = 16;  % # shards of A
opts.dir_scratch = dir_scratch;


% opts.stvec_flag = 1;   % 0 = regularization against rigid model (i.e.; starting value is not supplied by rc)
opts.distributed = 0;
opts.disableValidation = 1;
opts.edge_lambda = opts.lambda;
opts.use_peg = 0;

% % configure point-match filter
opts.verbose = 1;
opts.debug = 0;

% solve (mainly for translation here because the point matches are sparse)

% rcfine.stack = ['Test_slab_' num2str(nfirst) '_' num2str(nlast)...
%     '_fine_nbrs_' num2str(opts.nbrs) '_deg_1'];
% 
% rcfine.versionNotes = gen_versionNotes(opts);
% if stack_exists(rcfine), delete_renderer_stack(rcfine);end
% [err,R, Tout] = system_solve(nfirst, nlast, rcsource, pm, opts, rcfine);disp(err);



% obtain high-quality point-matches based on intensity-based registration

% read sections (tiles) in range and configure them
L = Msection(rc, nfirst, nlast);
for ix = 1:numel(L.tiles)
    L.tiles(ix).server = rc.baseURL;
    L.tiles(ix).project = rc.project;
    L.tiles(ix).stack   = rc.stack;
    L.tiles(ix).owner   = rc.owner;
    
    L.tiles(ix).fetch_local = 1;
end

% generate neighbor pair list
z_list = [L.tiles(:).z];
counter = 1;
pair_list = [];
for ix1 = 1:numel(z_list)   % assuming one tile per section
    z1_ix = find(z_list==z_list(ix1));  % index of tiles(section) with z value z_list(ix1);
    z     = z_list(ix1);
    for zix2 = z+1:z+opts.nbrs  % 
        z2_ix = find(z_list==zix2);
        if ~isempty(z2_ix)
            pair_list(counter,:) = [z1_ix z2_ix]; % indices into the L.tiles array
            counter = counter + 1;
        end
    end
end
%disp(pair_list);
tiles = L.tiles;

% process each pair to determine point-matches
kk_clock;
delete_point_match_collection(pmfine);  % we will start with a new point-match collection
disp(' ');
parfor_progress(size(pair_list,1));
parfor ix = 1:size(pair_list,1)
    % tile pair
    t1 = tiles(pair_list(ix,1));
    t2 = tiles(pair_list(ix,2));
    
    % get (render) images
    im1 = get_image(t1);
    im2 = get_image(t2);
    
    % reduce scale
    im1 = mat2gray(imresize(im1, pm_fine_scale));
    im2 = mat2gray(imresize(im2, pm_fine_scale));

    [tform, show_pair_image{ix},  p1, p2] = ...
        register_image_pair_intensity_based(im1, im2,...
        pmfine.npoints_per_pair, 0);
    p1 = p1./pm_fine_scale;
    p2 = p2./pm_fine_scale;

    % generate point-pair set struct, json and ingest
    ingest_point_match_set(pmfine, t1.renderer_id, t2.renderer_id, ...
                          t1.sectionId, t2.sectionId, p1, p2, dir_scratch);
    parfor_progress;
end
parfor_progress(0);
kk_clock;

%  show image registration residuals
error = [];

for ix = 1:size(pair_list,1)
    im = double(show_pair_image{ix});
    im(:,:,1) = mat2gray(im(:,:,1));%/norm(im(:,:,1));
    im(:,:,2) = mat2gray(im(:,:,2));%/norm(im(:,:,2));
    fh = figure('visible',visibility);
    clf;imshow(im);title([ num2str(pair_list(ix,:)) '----   Serial:' num2str(ix)]);drawnow;pause(1);
    saveas(fh, [alignmentDiagnosticsDir sprintf('alignment_%3.3d_%3.3d.tif', pair_list(ix,1), pair_list(ix,2))]);
    close(fh);
    error(ix) = sqrt(sum(sum((im(:,:,1)-im(:,:,2)).*(im(:,:,1)-im(:,:,2)))));
    disp(error(ix));
end
error = (error)/norm(error);
figure;plot(error, 'k*');
transopts = opts;

%  solve and ingest
if strcmp(fineOrFiner, 'fine')
    transopts.degree = 0;
    rctranslation = rcfine;
    rctranslation.stack = strrep(rctranslation.stack, 'fine','translation');
    rctranslation.versionNotes = 'Translation only';
    if stack_exists(rctranslation), delete_renderer_stack(rctranslation);end
    [err,R, Tout] = system_solve_translation(nfirst, nlast, rc, pmfine, transopts, rctranslation);
    rc = rctranslation;
end
rcfine.versionNotes = [ rcfine.versionNotes '. Parameters: ' gen_versionNotes(opts)];
if stack_exists(rcfine), delete_renderer_stack(rcfine);end
[err,R, Tout] = system_solve(nfirst, nlast, rc, pmfine, opts, rcfine);disp(err);
end

function generateSingleChannelScape( rc, root_directory, scale)

command = ['/groups/flyTEM/flyTEM/render/spark/render_scapes.sh ' ...
  '--sparkRunInBackground n --sparkNodes 1 --sparkBillTo loogerl --sparkDir /groups/looger/loogerlab/render/spark_output ' ... 
  '--baseDataUrl http://tem-services.int.janelia.org:8080/render-ws/v1 ' ...
  '--owner ' rc.owner ' --project ' rc.project ' --stack ' rc.stack ' ' ...
  '--rootDirectory ' root_directory ' --format png ' ...
  '--scale ' num2str(scale) ];
system(command); %can get output, and write output to spark_launch.log

end

function generateAllScapes( rc, channels, root_directory, scale)
% Loop through to generate all scapes
num_channels = numel(channels);
rc_stack_template = rc.stack;
for channel_index = 1:num_channels
    rc.stack = strrep(rc_stack_template, channels{1}, channels{channel_index});
    generateSingleChannelScape( rc, root_directory, scale);
end
fprintf('DONE!\n');
end

function createOutputImages( rc, channels, channelOrder, root_directory, scale, outputType)
% Generate the scapes and output figures
generateAllScapes( rc, channels, root_directory, scale);
channelOutputDirectories = cell(numel(channels),1);
for channelIndex=1:numel(channels)
    directories = dir(strrep(sprintf('%s%s/%s/', root_directory, rc.project, rc.stack), channels{1},channels{channelIndex}));
    channelOutputDirectories{channelIndex} = [directories(end).folder '/' directories(end).name '/000/'];
end
datetimeNow = datetime('now');
mergedOutputDirectory = strrep(channelOutputDirectories{1},channels{1},'merged');
mergedOutputDirectorySplit = strsplit(mergedOutputDirectory, '_');
mergedOutputDirectorySplit(end-1:end) = {datestr(datetimeNow, 'yyyymmdd'), datestr(datetimeNow, 'HHMMSS')};
mergedOutputDirectory = [strjoin(mergedOutputDirectorySplit, '_') '/000/'];
system(['mkdir -p ' mergedOutputDirectory]);
files = dir([channelOutputDirectories{1} '*.png']);
parfor currentFileIndex = 1:numel(files)
    outputImage = [];
    for channelIndex=1:numel(channels)
        currentChannelFileName = [channelOutputDirectories{channelIndex} files(currentFileIndex).name];
        outputImage(:,:,channelOrder(channelIndex)) = imread(currentChannelFileName);
    end
    if strcmp(outputType,'double')
        imwrite(double(outputImage), [mergedOutputDirectory strrep(files(currentFileIndex).name,'.png', '.tif')]);%,'Compression', 'none');
    elseif strcmp(outputType,'uint8')
        imwrite(uint8(outputImage), [mergedOutputDirectory strrep(files(currentFileIndex).name, '.png', '.tif')]);%,'Compression', 'none');
    elseif strcmp(outputType, 'uint16')
        imwrite(uint8(outputImage), [mergedOutputDirectory strrep(files(currentFileIndex).name, '.png', '.tif')]);%,'Compression', 'none');
    else
        error(['Unrecognized output type ' outputType])
    end
end
end