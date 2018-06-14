function [Result] = diamond_vanish(InputImg, Normalization, SpaceSize, PatchSize, VanishNumber)
% Detection of vanishing point using diamond space.
%
% Results = DIAMOND_VANISH(InputImg, Normalization, SpaceSize, PatchSize)
%
% Input:
%  InputImg       input image
%  Normalization  float number; 
%                 normalization of the image (1 means normalization form -1 to 1)
%  SpaceSize      int number; 
%                 resolution of the accumulation space (final space has dims SpaceSize x SpaceSize)
%  PatchSize      int row vector; 
%                 radius of a patch, from which edge pixels are extracted and use for ellipse fitting for detection of orientation of the edge point
%
% Output:
%  Results structure with fields:
%    Results.Space         accumulated diamond space (further is used for orthogonalization)
%    Results.PC_VanP       positions of the maxima in R.Space
%    Results.PC_VanP_Norm  normalized position of the maxima (R.Space bounded from -1 to 1)
%    Results.CC_VanP       position of the vanishing point in the input image coordinates


%find edge points
EdgeImg = edge(rgb2gray(InputImg),'canny'); 

%find lines
LinesData = mx_lines(int32(padarray(EdgeImg,[PatchSize(end),PatchSize(end)])), int32(PatchSize)); 
LinesData(3,:) = LinesData(3,:)*Normalization;

SubPixelRadius = 2;
Threshold = 0.05;
Result = struct('Space', [], 'PC_VanP',[], 'PC_VanP_Norm', [], 'CC_VanP', []);

%diamond stuff
for V = 1:VanishNumber           
    %raster and find space
    space = mx_raster_space(SpaceSize, LinesData);
    %figure
    %imagesc(space);
    Result.PC_VanP(V,:) = find_maximum(space, SubPixelRadius);
    if(V == 1) Result.Space = space; end
          
    %normalize result
    Result.PC_VanP_Norm(V,:) = normalize_PC_points(Result.PC_VanP(V,:), SpaceSize);
        
    %get lines close to VP
    Distance = point_to_lines_dist(Result.PC_VanP_Norm(V,:), LinesData(1:3,:)');
        
    %remove lines 
    LinesData(:,(Distance < Threshold)') = [];        
end
    
Result.CC_VanP = PC_point_to_CC(Normalization, Result.PC_VanP_Norm, size(InputImg));

end

function NormVP = normalize_PC_points(VanP, SpaceSize)
    NormVP = (2.*VanP -(SpaceSize + 1))./(SpaceSize - 1);
end

function D = point_to_lines_dist(Point, Lines)
    x = Point(1);
    y = Point(2);

    T = [0,-1,1;1,-1,0;0,-1,0;...
         0,-1,1;1,1,0 ;0,-1,0;...
         0,1,1 ;1,-1,0;0,-1,0;...
         0,1,1 ;1,1,0 ;0,-1,0];
    
    L = Lines*T';
    
    P(:,1) = (L(:,1:3)*[x,y,1]')./sqrt(sum(L(:,1:2).^2,2));
    P(:,2) = (L(:,4:6)*[x,y,1]')./sqrt(sum(L(:,4:5).^2,2));
    P(:,3) = (L(:,7:9)*[x,y,1]')./sqrt(sum(L(:,7:8).^2,2));
    P(:,4) = (L(:,10:12)*[x,y,1]')./sqrt(sum(L(:,10:11).^2,2));
   
    D = min(abs(P),[],2);
end


function [VanPC] = find_maximum(Space, R)
    [r, c] = find(max(Space(:)) == Space);

    S = padarray(double(Space),[R,R]);    
    
    O = S(r(1):r(1)+R*2, c(1):c(1)+R*2);
    
    [mc,mr] = meshgrid([-R:R],[-R:R]);
    SR = O.*mr;
    SC = O.*mc;
    
    C = c(1) + sum(SC(:))/sum(O(:));
    R = r(1) + sum(SR(:))/sum(O(:));
    
    VanPC = [C,R];
end

function VanishCC = PC_point_to_CC(Normalization, VanishPC, ImgSize)
    u = VanishPC(:,1);
    v = VanishPC(:,2);
    NormVanishCC = [v, sign(v).*v + sign(u).*u - 1, u];
    
    reg = abs(NormVanishCC(:,3)) > 0.005;
    
    NormVanishCC(reg,:) = bsxfun(@rdivide, NormVanishCC(reg,:), NormVanishCC(reg,3));
    NormVanishCC(~reg,:) = normr(NormVanishCC(~reg,:));
    
    VanishCC = NormVanishCC;
    VanishCC(~reg,3) = 0;
    
    w = ImgSize(2);
    h = ImgSize(1);
    m = max(ImgSize);
    
    VanishCC(reg,1) = (VanishCC(reg,1)./Normalization.*(m - 1) + w + 1)./2;
    VanishCC(reg,2) = (VanishCC(reg,2)./Normalization.*(m - 1) + h + 1)./2;
end
