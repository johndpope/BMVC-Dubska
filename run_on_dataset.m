function Results = run_on_dataset(DB, DB_path, Normalization, SpaceSize, PatchSize)
% Detection of vanishing point using diamond space and computation of world directions (both nonorthogonal and orthogonal)
%
% Results = RUN_ON_DATASET(DB, DB_path)
% Results = RUN_ON_DATASET(DB, DB_path, Normalization, SpaceSize, PatchSize)
%
% Input:
%  DB             cell array with images names (e.g. Manhattan_Image_DB_Names)
%  DB_path        path to dataset data (with camera parameters and image data)%
%  Normalization  float number; 
%                 normalization of the image (1 means normalization form -1 to 1)
%  SpaceSize      int number; 
%                 resolution of the accumulation space (final space has dims SpaceSize x SpaceSize)
%  PatchSize      int row vector; 
%                 radius of a patch, from which edge pixels are extracted and use for ellipse fitting for detection of orientation of the edge point
%
%  default parameters are from training on a ECCV_TrainingAndTestImageNumbers -> trainingSetIndex	
%    Normalization = 0.4
%    SpaceSize = 321
%    PatchSize = [6:4:22]
%
% Output:
%  Results structure with fields:
%    Results.Image        image name
%    Results.Detected     3x3 matrix with three detected vanishing points, each in one column
%    Results.Orthogonal   3x3 matrix with three vanishing points after orthogonalization, each in one column

    Camera = load([DB_path '/cameraParameters.mat']);

    if ~exist('Normalization','var') || isempty(Normalization) Normalization = 0.4; end
    if ~exist('SpaceSize','var') || isempty(SpaceSize) SpaceSize = 321; end
    if ~exist('PatchSize','var') || isempty(PatchSize) PatchSize = [6:4:22]; end
    
    for i = 1:length(DB)
        img = DB{i};
        
        %load image and GT
        I = imread([ DB_path '/' img(1:end-1) '/' img(1:end-1) '.jpg']);
        
        %detect vanish
        Detect = diamond_vanish(I, Normalization, SpaceSize, PatchSize, 3);
        
        %find ortho
        Results(i).Image = img;
        Results(i).Detected = img_point_to_world(Detect.CC_VanP, Camera.focal, Camera.pp, Camera.pixelSize)';
        Results(i).Orthogonal = find_ortho(Results(i).Detected', Camera, Normalization, SpaceSize, Detect.Space, size(I))';
    end
    
end


function OV = find_ortho(DetectedVP, Camera, Normalization, SpaceSize, Space, ImgSize)
    V = zeros(9,3);    
    
%     [V(1:3,:), V(4:6,:), V(7:9,:)] = make_all_ortho(DetectedVP);
    V(1:3,:) = make_ortho(DetectedVP);
    V(4:6,:) = make_ortho(DetectedVP([2,1,3],:));
    V(7:9,:) = make_ortho(DetectedVP([3,1,2],:));
    
    ImgPoints = world_point_to_img(V, Camera.focal, Camera.pp, Camera.pixelSize);
    
    w = ImgSize(2);
    h = ImgSize(1);
    m = max([w,h]);
        
    reg = abs(ImgPoints(:,3)) > eps;
        
    NormOneR = ImgPoints(:,2);
    NormOneC = ImgPoints(:,1);
        
    NormOneR(reg) = (2*NormOneR(reg) - (h + 1))/(m - 1) .* Normalization;
    NormOneC(reg) = (2*NormOneC(reg) - (w + 1))/(m - 1) .* Normalization;
    
    Ortho_Detect = CC_point_to_PC(SpaceSize, [NormOneC, NormOneR, ImgPoints(:,3)]);
    Vals = interp2(single(Space), Ortho_Detect(:,1), Ortho_Detect(:,2));
    Vals([1,4,7]) = Vals([1,4,7])./Vals(1);
    Vals([2,5,8]) = Vals([2,5,8])./Vals(5);
    Vals([3,6,9]) = Vals([3,6,9])./Vals(9);
              
    [~,I] = max([sum(Vals(1:3)),sum(Vals(4:6)),sum(Vals(7:9))]);
    OV =  V(I*3 - 2: I*3,:); 
    
end

function ImgPoint = world_point_to_img(WorldPoint, Focal, PrincipalPoint, PixelSize)    
    reg = abs(WorldPoint(:,3)) > 0.005;
    
    WorldPoint(reg,:) = bsxfun(@rdivide, WorldPoint(reg,:), WorldPoint(reg,3));
    WorldPoint(~reg,3) = 0;
    
    k = PixelSize./Focal;
    
    ImgPoint = [WorldPoint(:,1)./k + PrincipalPoint(1)*WorldPoint(:,3), -WorldPoint(:,2)./k + PrincipalPoint(2)*WorldPoint(:,3), WorldPoint(:,3)];        
end


function WorldPoint = img_point_to_world(ImgPoint, Focal, PrincipalPoint, PixelSize)
    k = PixelSize./Focal;
    WorldPoint = [(ImgPoint(:,1) - PrincipalPoint(1)*ImgPoint(:,3)).*k, -(ImgPoint(:,2) - PrincipalPoint(2)*ImgPoint(:,3)).*k,ImgPoint(:,3)];
    WorldPoint = normr(WorldPoint);
end

function V = make_ortho(Van) 
	V(1,:) = Van(1,:);    
    V(2,:) = Van(2,:) - (Van(2,:)*Van(1,:)')*Van(1,:);
	V(3,:) = Van(3,:) - (Van(3,:)*Van(1,:)')*Van(1,:);
   
    if(abs(V(2,:)*Van(2,:)') > abs(V(3,:)*Van(3,:)'))
        V(3,:) = cross(V(1,:),V(2,:)); 
    else
        V(2,:) = cross(V(1,:),V(3,:)); 
    end            
end

function [VanishPC, VanishPC_S] = CC_point_to_PC(SpaceSize, VanishCC)
    x = VanishCC(:,1);
    y = VanishCC(:,2);
    w = VanishCC(:,3);
    NormVanishPC = [-w, -x, sign(x.*y).*x + y + sign(y).*w];
    
    VanishPC_S = bsxfun(@rdivide, NormVanishPC, NormVanishPC(:,3));
    VanishPC = (VanishPC_S(:,1:2) .* (SpaceSize - 1) + (SpaceSize + 1))./2;
end


