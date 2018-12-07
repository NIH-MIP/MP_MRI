%% Align two sets of dicom images based on dicom coordinates
%    performs affine alignment and linear interpolation to align series
%   
%   this stand alone version reads in full DICOM folder for an MRI exam and
%   automatically aligns high b-value and ADC maps to T2 imaging
%
%   input:
%       1. parent DICOM folder for MRI exam
%
%   output: 
%       1. aligned highB 
%               original untransformed saved to ./highB/raw/
%               new transformed saved to ./highB/aligned/
%       2. aligned ADC 
%               original untransformed saved to ./ADC/raw/
%               new transformed saved to ./ADC/aligned/ 
%
% Created by Stephanie Harmon 10/2018
% Portions based on work by Daniel R Warren http://github.com/drw25                                    
%
% Requires MATLAB Image Processing toolbox.


function [] = seriesAlign3D(varargin)

    %parse through input or select directory
        if length(varargin) == 1 
            directory = varargin{1};   
        else 
            disp('You must provide a parent DICOM directory');
            return
        end

    %check to see if T2 (reference) folder exists, if not -> exit
        if(~exist([directory '\t2']))
            disp('WARNING! No T2 folder. Exiting');
            return
        else
            reference_loc = [directory '\t2'];
        end

    %begin alignment - identify ADC series     
        if(~exist([directory '\adc']))
            disp('WARNING! No ADC folder'); %send warning if not there
        else
            %if its there - send to runAlignment
            target_loc = [directory '\adc'];
            disp(['transforming ADC']);
            runAlignment(reference_loc, target_loc);
        end        
        
    %begin alignment - identify highB series    
        if(~exist([directory '\highb']))
            disp('WARNING! No highB folder'); %send warning if not there
        else
            %if its there - send to runAlignment 
            target_loc = [directory '\highb'];
            disp(['transforming high B']);
            runAlignment(reference_loc, target_loc);
        end
        

end




%% runAlignment
function [] = runAlignment(reference_loc, target_loc)
% this function runs the alignment after calling DICOM transformation
% and saves to new folders following writing new DICOMS after alignment
% input: DICOM folders of reference and target DICOMS
% output: SPECIFIC TO DICOM DATABASE FOR MIP
%           creates new folder "raw" for old DICOMS
%           saves aligned images as "aligned" 
%           reference folder remains unchanged
    %grab coordinate data and parse dicoms
        %disp('loading target series to be transformed... ');
        [r_img r_ref r_trnsf r_dcm] = ReadDicom3D(reference_loc, 'F');
        %disp('loading reference series for transformation... ');
        [t_img t_ref t_trnsf t_dcm] = ReadDicom3D(target_loc, 'T');

    % affine alignment of target to reference
        r_img = double(r_img);
        t_img = double(t_img);
        t_to_r = t_trnsf/r_trnsf;
        t_to_r(:,4) = [0;0;0;1]; % fix precision errors
        t_to_r_trsnf = affine3d(t_to_r);
        t_aligned = imwarp(t_img,t_ref,t_to_r_trsnf,'OutputView',r_ref);

    %save new dicom
        saveFolder = [target_loc '\aligned'];
        mkdir(saveFolder);
        %disp('saving aligned images... ');
        matchdicoms(r_dcm, t_dcm, t_aligned, saveFolder)

end





%% readDicom3D
function [img ref_data trsf_data dcm_full] = ReadDicom3D(folder, move)
% this function reads in dicom info and outputs various forms of reference
% data and transformation data
%
% INPUTS:   folder containing DICOM images
%
%
% OUTPUTS:
%           img  : image data (matrix)
%           ref_img: a frame of reference in IMAGE SPACE but with PHYSICAL
%                  UNITS (ie. millimetres not pixels) in imref3d format
%           T_img: a 4D transformation matrix suitable for use with imwarp
%                  that maps image space to the DICOM reference coordinate
%                  system
%           dmeta: Subset of DICOM metadata. UIDs and 3 important geometric
%                  fields are retained for each slice (IOP, IPP, PS).
%                  Remainder is only first-level fields whose value is
%                  identical for all instances in the series
%
% NOTE THIS CODE IS FOR 3D DATA ONLY - ASSUMES ALL HAVE SAME Z
%           CURRENTLY 
% 
% PARSE IMAGE DIRECTORY
    listing = dir(folder); 
    num_slices = 0;
    zpos = zeros(1,1);
    switch move
        case 'T'
            mkdir([folder '\raw']);
        case 'F'
    end
        
    
    %iterate each slice and grab data
    for j = 1:numel(listing)
        file_j = [folder filesep listing(j).name];
        if ~listing(j).isdir && isdicom(file_j)
            num_slices = num_slices+1;
            dinfo{num_slices} = dicominfo(file_j); % grab all dicom metadata
            zpos(num_slices) = dinfo{num_slices}.ImagePositionPatient(3);
            img(:,:,num_slices) = dicomread(dinfo{num_slices}); % read in slice
            switch move
                case 'T'
                    movefile(file_j,[folder '\raw']);
                case 'F'
            end
        end
    end
    
    %sort images and dicom in ascending order
    [B,ordering] = sort(zpos,2); %B is unused dummy variable
    img = img(:,:,ordering);
    dinfo = dinfo(ordering);
    
    %get PATIENT position info - grab from top slice, will be same for all slices 
        IOP_xyz = dinfo{1}.ImageOrientationPatient; % direction cosines relating i vector to xyz and j vector to xyz
        IPP_xyz = dinfo{1}.ImagePositionPatient; % (x,y,z) position of (i,j)=(0,0)
        idir_xyz = IOP_xyz(1:3); % direction vector for rows in xyz
        jdir_xyz = IOP_xyz(4:6); % direction vector for columns in xyz
        kdir_xyz = cross(idir_xyz,jdir_xyz); % previously compared IPP(1) and IPP(end) but sign problems in sagittal images
        %normalize by vector magnitude (note to self, usually = 1)
            idir_xyz = idir_xyz/norm(idir_xyz); 
            jdir_xyz = jdir_xyz/norm(jdir_xyz);
            kdir_xyz = kdir_xyz/norm(kdir_xyz);
            
    %get SLICE (k) and PLANE (i,j) info for all voxels [note plane same in all slices)  
        %get vector of SLICE position info 
            img_vk = cellfun(@(x)x.ImagePositionPatient,dinfo,'UniformOutput',false); % img_vk will be vector of CENTRAL voxel in slice k
            img_vk = bsxfun(@minus,[img_vk{:}],dinfo{1}.ImagePositionPatient); %re-orient all to first slice
            img_vk = sqrt(sum(img_vk.^2,1)); %notice here we now have the z oriented with correct distance between planes (thickness)
        %get spacing info for 3D volume        
            img_di = dinfo{1}.PixelSpacing(1); % row resolution (voxel size)
            img_dj = dinfo{1}.PixelSpacing(2); % column resolution (voxel size)
            img_dk = img_vk(2)-img_vk(1); %slice thickness
                %note to self, we already knew this but now its derived from data itself 
                %instead of relying on slice thickness and forcing data backwards to fit
        %get vector of voxel position within each slice   
            img_vi = img_di*double(0:(-1+dinfo{1}.Rows)); % img_vi is vector of pixel centres in i
            img_vj = img_dj*double(0:(-1+dinfo{1}.Columns)); % img_vj is vector of pixel centres in j

    
    % ROTATE TO PATIENT COORDINATES AND TRANSLATE TO WORLD COORDINATES
        rotmat = [idir_xyz jdir_xyz kdir_xyz]; % Rotation matrix to map ijk directions to patient xyz directions
        rotmat(4,4) = 1;
        tmat = [[eye(3); IPP_xyz'] [0;0;0;1]]; % Translation matrix to map origin to IPP(0)
        
    % physical FOV location
        img_fovi = [img_vi(1)-img_di/2 img_vi(end)+img_di/2];
        img_fovj = [img_vj(1)-img_dj/2 img_vj(end)+img_dj/2];
        img_fovk = [img_vk(1)-img_dk/2 img_vk(end)+img_dk/2];

    % IMREF3D is a beautiful thing. look it up.
    % puts dicom into frame of reference based on physical FOV location
        ref_data = imref3d(size(img),img_fovj,img_fovi,img_fovk); % Frame of reference for image space
        trsf_data = rotmat'*tmat; % Transformation matrix from image coordinates to real space
             % Matrix is right-multiplied by imwarp, so earlier operations should appear first

    %SEND BACK ORDERED INFO FOR EASIER REWRITING OF DICOMS OUT        
        dcm_full = dinfo;

        
%     %ALSO SEND BACK META FOR WRITING TO NEW FILES
%         dcm_meta = struct();
% 
%         % Keep a subset of DICOM metadata - only first-level fields with numeric or
%         % char datatypes, whose value is the same for every instance
%         % note from SH: didnt change anything in below code from original
%         file_j = fieldnames(dinfo{1});
%         for i = 1:numel(file_j)
%             keep = false;
%             refinfo = dinfo{1}.(file_j{i});
%             for j = 2:numel(dinfo)
%                 if isfield(dinfo{j},file_j{i})
%                     testinfo = dinfo{j}.(file_j{i});
%                     if strcmp(class(refinfo),class(testinfo)) && ... 
%                         ( ...
%                          ( isnumeric(refinfo) && (numel(refinfo) == numel(testinfo)) && all(refinfo == testinfo) ) ...
%                         || ...
%                          ( ischar(refinfo) && strcmp(refinfo,testinfo) ) ...
%                         )
%                             keep = true;
%                     end
%                 else
%                     keep = false;
%                 end
%                 if ~keep; break; end;
%             end
%             if keep
%                 dcm_meta.(file_j{i}) = dinfo{1}.(file_j{i});
%             end
%         end
% 
%         % Keep some important DICOM metadata for every slice
% 
%         fieldcontents = @(f)cellfun(@(x){x.(f)},dinfo);
% 
%         % SH note: this was from old code, I scrapped it
%         %
%         % dicom_meta.SOPInstanceUID = fieldcontents('SOPInstanceUID');
%         % dicom_meta.ImageOrientationPatient =  fieldcontents('ImageOrientationPatient');
%         % dicom_meta.ImagePositionPatient = fieldcontents('ImagePositionPatient');
%         % dicom_meta.PixelSpacing = fieldcontents('PixelSpacing');

end





%% matchdicoms
function [dcm_matched] = matchdicoms(dcm_ref, dcm_orig, img_align, saveFolder)
% Created by Stephanie Harmon 9-28-2018
% 
% Write aligned image series to new dicom data
%    This code is intended to be used after alignment to world coordinates 
%    and resampling to reference size (i.e. PET --> CT or ADC --> T2W)
%
% INPUTS:
%   1. dcm_ref = array of dicom structure data from reference sequence 
%                (i.e. CT or T2W) from full image series 
%
%   2. dcm_orig = single file meta data captured from ReadDicom3D of
%                 aligned image set (i.e. PET or ADC in this example). this
%                 will be incomplete for new data (patient orientation,
%                 image position, and slice resolution/thickness)
%
%   3. img_align = matrix of aligned image data to be exported to dicom
%
%   4. saveFolder = directory to save new dicom data of aligned set
    %mkdir(saveFolder); 

    %global variables (i.e. those that dont change by slice)
    %IOP, resolution, and size information

    
    %slice-specific variables 
    %note in future versions i should probably call back an error if
    %dcm_ref array size is different from the number of slices in img_align
        for k = 1:size(img_align,3)
            if(k>length(dcm_orig))
                dcm_meta = dcm_orig{length(dcm_orig)};
                MS_UID = dcm_meta.MediaStorageSOPInstanceUID;
                MS_UID_digit = str2num(MS_UID(end))+1;
                MS_UID(end) = MS_UID_digit;
                dcm_meta.MediaStorageSOPInstanceUID = MS_UID;
                dcm_meta.SOPInstanceUID = MS_UID;
            else
                dcm_meta = dcm_orig{k};
            end
            dcm_meta.ImageOrientationPatient = dcm_ref{1}.ImageOrientationPatient;
            dcm_meta.PixelSpacing = dcm_ref{1}.PixelSpacing;
            dcm_meta.Width = dcm_ref{1}.Width;
            dcm_meta.Height = dcm_ref{1}.Height;
            dcm_meta.Rows = dcm_ref{1}.Rows;
            dcm_meta.Columns = dcm_ref{1}.Columns;
            dcm_meta.PixelSpacing = dcm_ref{1}.PixelSpacing;
            dcm_meta.SliceThickness = dcm_ref{1}.SliceThickness;
            dcm_meta.SliceLocation = dcm_ref{k}.SliceLocation;
            sliceData = uint16(img_align(:,:,k)); %convert back to unsigned int
            dcm_meta.ImagePositionPatient = dcm_ref{k}.ImagePositionPatient; %IPP 
            dicomwrite(sliceData,[saveFolder '\slice' int2str(k) '.dcm'],dcm_meta); %save
        end
    
end