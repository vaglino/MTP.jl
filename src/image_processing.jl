#= 
image_processing.jl contains the main functions to process TIRF microscope images, segment and extract single cells
Stefano Travaglino, 2023, Zhu Lab, Georgia Tech
=#

# module ImageProcessing


# using TiffImages, Images, Statistics, TiffImages, Revise
# includet("background.jl")
# using .IlluminationCorrection
# using PyCall
# morph = pyimport("skimage.morphology")
# feat = pyimport("skimage.feature")
# flt = pyimport("skimage.filters")
# ndimage = pyimport("scipy.ndimage")
# skseg = pyimport("skimage.segmentation")
# skexposure = pyimport("skimage.exposure")



# export Cell, files2cells, img2cells

mutable struct Cell
    I
    seg
    t
    coord
    I_fixed
    seg_fixed
end

"""
Main function to load and process files into single cells
# Arguments
- `data_dir`: directory containing tif files
- `tifs_names`: names of tif files to be processed
- `illumination`: average illumination of image stack
# Returns
- `cells_stack`: stack of cells
"""
function files2cells(data_dir,tifs_names;illumination=nothing)

    illumination_smooth = imfilter(illumination, Kernel.gaussian(15))
    cells_stack = []
    for tif_name in tifs_names
        @show tif_name
        img = Float64.(TiffImages.load(data_dir*tif_name)) 

        # uncomment if only two channels
        # img_new = cat(rand(Float64,size(img)[1:2]) , img, dims=3)
        # img_new = img_new[:,:,[1,2,3]]

        # @show size(img_new)

        # uncoment if three channels
        img_new = img[:,:,[1,2,3]]


        @time cells, I_sub, DIC_labels, overlays  = img2cells(img_new,illumination,illumination_smooth)

        push!(cells_stack, cells)
    end
    return cells_stack
end

"""
From a single image, extract cells
# Arguments
- `img`: image to be processed
# Returns
- `cells`: array of cells
"""
function img2cells(img)
    plotting = false
    _I_sub, gauss = subtract_background(img) .* 2^16 # CHECK IF IMAGE IS 16 bit or 12 BIT!!!!!
    # imshow(I_sub[3])
    I_sub = arr2tens(_I_sub)
    # imshow(I_sub[:,:,3])
    I_sub[:,:,2] = img[:,:,2] ./ illumination_RICM

    I_DIC_rsc = skexposure.rescale_intensity(I_sub[:,:,3]) 
    # imshow(I_DIC_rsc)

    I_seg_DIC, I_bw_DIC = segment_DIC(I_DIC_rsc)
    # I_seg_DIC_color = plot_segmented(I_seg_DIC, I_bw_DIC)
    # imshow(I_seg_DIC_color)

    DIC_labels = labels_map(I_seg_DIC) .* I_bw_DIC
    # imshow(DIC_labels)

    cells = seg2cell(I_sub,DIC_labels)  
    overlay = plot_segmentation_overlay(I_DIC_rsc, DIC_labels)

    plotting && imshow(overlay)

    _cells = add_RICM_segmentation_to_cell.(cells)
    overlays = mosaic_cropped_cells(_cells)

    return _cells, I_sub, DIC_labels, overlays  
end

function img2cells(img,illumination,illumination_smooth;MTP_channel=1)
    plotting = true

    I_MTP = img[:,:,MTP_channel] .* 2^16
    # I_MTP_no_outl = remove_outliers(I_MTP;radius=15,threshold=2.5)
    # I_MTP = illumination_correction(I_MTP_no_outl,illumination)
    I_MTP = illumination_correction_MTP(I_MTP,illumination,illumination_smooth)

    # imshow(I_MTP)
    # _I_sub, gauss = subtract_background(img;channels=[2,3,4]) .* 2^16 
    # _I_sub, gauss = subtract_background(img;channels=[2,3]) .* 2^16 
    # imshow(I_sub[3])
    # I = push!(_I_sub,I_MTP)
    # I = I[[3,1,2]]
    # I = I[[3,1,2]]
    I = [I_MTP, img[:,:,2] .* 2^16, img[:,:,3] .* 2^16]


    I_sub = arr2tens(I)
    # imshow(I_sub[:,:,3])
    I_sub[:,:,2] = img[:,:,2] ./ illumination_RICM
    I_sub[:,:,2] = skexposure.rescale_intensity(I_sub[:,:,2])
    I_sub[:,:,2] = skexposure.equalize_adapthist(I_sub[:,:,2])


    I_DIC_rsc = skexposure.rescale_intensity(I_sub[:,:,3]) 
    # imshow(I_DIC_rsc)

    I_seg_DIC, I_bw_DIC = segment_DIC(I_DIC_rsc)
    # I_seg_DIC_color = plot_segmented(I_seg_DIC, I_bw_DIC)
    # imshow(I_seg_DIC_color)

    DIC_labels = labels_map(I_seg_DIC) .* I_bw_DIC
    # imshow(DIC_labels)

    cells = seg2cell(I_sub,DIC_labels)  
    overlay = plot_segmentation_overlay(I_DIC_rsc, DIC_labels)
    plotting && imshow(overlay)

    _cells = add_RICM_segmentation_to_cell.(cells)
    overlays = mosaic_cropped_cells(_cells)

    return _cells, I_sub, DIC_labels, overlays  
end


function seg2cell(I,seg)
    # for every segmented cell in an image output all cells structs
       
    # I_ind =[ [i,j] for i in 1:size(I,1), j in 1:size(I,2)]
    cells = []
    for i in 1:maximum(seg)

        # crop a single cell
        I_cell,segmented = crop_cell(I,seg,i)
        # repeat cropping, with fixed size, for representative Images
        I_cell_fixed,segmented_fixed = crop_cell_fixed_size(I,seg,i)

        t = 0
        coord_cell = [] #I_ind[ lb_y-δ : ub_y+δ , lb_x-δ : ub_x+δ ] 
        # push!(cells, Cell1(I_cell,segmented,t,coord_cell))
        push!(cells, Cell(I_cell,segmented,t,coord_cell,I_cell_fixed,segmented_fixed))
    end

    return cells
end

 



"""
This function takes a Cell object and adds the RICM segmentation to it
# Arguments
- `cell`: Cell object
# Returns
- `cell`: Cell object with RICM segmentation added
"""
function add_RICM_segmentation_to_cell(cell)
    I_new = cell.I[:,:,2]
    seg, bw, all_segs = segment_RICM_single_cell(I_new)
    # imshow(bw)
    # @show seg
    # RICM_labels = labels_map(seg) .* bw
    # # imshow(RICM_labels)
    # overlap = RICM_labels .* (cell.seg .!= 0)

    # imshow(overlap)
    RICM_seg = Int.(bw)
    cell.seg = cat(cell.seg, RICM_seg, dims=3)
    
    MTP = cell.I[:,:,1]
    # @show sum(isnan.(MTP))
    backg = MTP .* (.!all_segs .| .!bw)
    backg_ = backg[backg .!= 0.0]
    # imshow(backg_)
    backg_val = mean(backg_) #+ 1*std(backg_)
    # @show mean(MTP)
    # @show backg_val
    MTP_sub = MTP .- backg_val
    # @show mean(MTP_sub)
    MTP_sub[MTP_sub .< 0] .= 0.0
    cell.I[:,:,1] = MTP_sub

    MTP_fixed = cell.I_fixed[:,:,1]
    MTP_fixed_sub = MTP_fixed .- backg_val
    MTP_fixed_sub[MTP_fixed_sub .< 0] .= 0.0
    cell.I_fixed[:,:,1] = MTP_fixed_sub

    return cell
end

"""
Plotting of the RICM segmentation results for all single cells
# Arguments
- 'cells': vector of Cell objects
# Returns
- 'overlays': vector of RICM segmentation overlay plots
"""
function mosaic_cropped_cells(cells;channel = 2)
    overlays = []
    for cell in cells
        I_cell = cell.I[:,:,channel]
        labels = cell.seg[:,:,2]
        # labels = morph.binary_erosion(labels,morph.disk(5))
        overlay = plot_segmentation_overlay(I_cell, labels)
        push!(overlays, overlay)
    end
    overlays
end

# end

using StaticArrays



function crop_cell(I, seg, i; crop_ratio=1.5)
    # Find all indices for the segmented cell i
    indxs = findall(x -> x == i, seg)
    if isempty(indxs)
        error("No pixels found for segment $i")
    end

    # Extract boundaries of the segmented cell
    indx_mat = hcat(getindex.(indxs, 1), getindex.(indxs, 2))
    ub_y = maximum(indx_mat[:, 1])
    lb_y = minimum(indx_mat[:, 1])
    ub_x = maximum(indx_mat[:, 2])
    lb_x = minimum(indx_mat[:, 2])

    # Calculate dimensions of the segmented cell
    y_dim = ub_y - lb_y + 1
    x_dim = ub_x - lb_x + 1

    # Determine expansion amounts based on the crop ratio
    expansion_y = round(Int, (y_dim / 2) * crop_ratio)
    expansion_x = round(Int, (x_dim / 2) * crop_ratio)

    # Compute desired expanded boundaries
    ub_Y = ub_y + expansion_y
    lb_Y = lb_y - expansion_y
    ub_X = ub_x + expansion_x
    lb_X = lb_x - expansion_x

    # Clip boundaries to stay within image dimensions
    _ub_Y = min(ub_Y, size(I, 1))
    _lb_Y = max(lb_Y, 1)
    _ub_X = min(ub_X, size(I, 2))
    _lb_X = max(lb_X, 1)

    # Initialize arrays for the cropped cell and segmentation
    desired_size_y = ub_Y - lb_Y + 1
    desired_size_x = ub_X - lb_X + 1
    I_cell = zeros(eltype(I), desired_size_y, desired_size_x, size(I, 3))
    seg_cell = zeros(eltype(seg), desired_size_y, desired_size_x)

    # Calculate offsets to place the cropped region within the initialized arrays
    offset_y = _lb_Y - lb_Y
    offset_x = _lb_X - lb_X

    # Crop the relevant portion from the image and segmentation
    cropped_I = I[_lb_Y:_ub_Y, _lb_X:_ub_X, :]
    cropped_seg = seg[_lb_Y:_ub_Y, _lb_X:_ub_X]

    # Place the cropped regions into the initialized arrays
    I_cell[offset_y+1:offset_y+size(cropped_I, 1), offset_x+1:offset_x+size(cropped_I, 2), :] = cropped_I
    seg_cell[offset_y+1:offset_y+size(cropped_seg, 1), offset_x+1:offset_x+size(cropped_seg, 2)] = cropped_seg

    # Zero out segmentation not equal to cell i
    seg_cell[seg_cell .!= i] .= 0

    return I_cell, seg_cell
end

"""
Crop a single cell with fixed size
# Arguments
- `I`: image to be processed
- `seg`: segmentation of image
- `i`: index of cell to be cropped
# Returns
- `I_cell`: cropped cell
"""
function crop_cell_fixed_size(I,seg,i)

    # crop_ratio = 2 # crop length in each direction compared to segmented
    # crop_size = round(6*6/√2 / px_size)  # ~25μm / px_size
    # crop_size = round(25 / px_size)  
    crop_size = round(Int, 17.5 / px_size)

    I_ind =[ [m,n] for m in 1:size(I,1), n in 1:size(I,2)]
    indxs = findall(x->x == i,seg)

    indx_mat = hcat(getindex.(indxs, 1), getindex.(indxs,2))

    # boundaries of segmented cell
    ub_y = maximum(indx_mat[:,1]) 
    lb_y = minimum(indx_mat[:,1])
    ub_x = maximum(indx_mat[:,2])
    lb_x = minimum(indx_mat[:,2])

    # dimensions of segmented cell
    y_dim = (ub_y-lb_y) + 1
    x_dim =  (ub_x-lb_x) + 1

    # δ is half the difference between predefined crop size and the cell size
    δ_y = Int(round( (crop_size-y_dim) /2))
    δ_x = Int(round( (crop_size-x_dim) /2))
    # if cell is larger than predefined crop size, expand by δ=1
    (δ_y > 0) || (δ_y = 1)
    (δ_x > 0) || (δ_x = 1)

    # initialize cropped cell (including zero-padded area)
    I_cell = zeros( y_dim+δ_y*2 , x_dim+δ_x*2 , size(I,3))
    seg_cell = zeros( y_dim+δ_y*2 , x_dim+δ_x*2 )

    # expand boundaries by δ
    ub_Y = ub_y + δ_y
    lb_Y = lb_y - δ_y
    ub_X = ub_x + δ_x
    lb_X = lb_x - δ_x

    # use smallest bounds between expanded bounds and size limit of image
    _ub_Y = min(ub_Y, size(I,1) )
    _lb_Y = max(lb_Y, 1)
    _ub_X = min(ub_X, size(I,2) )
    _lb_X = max(lb_X, 1)

    # crop image 
    cropped_I = I[_lb_Y:_ub_Y , _lb_X:_ub_X , :]
    cropped_seg = seg[_lb_Y:_ub_Y , _lb_X:_ub_X]

    # assign cropped cell to initialized cell (which includes zero-padded area)
    I_cell[_lb_Y-lb_Y+1 : -lb_Y+_ub_Y+1,
           _lb_X-lb_X+1 : -lb_X+_ub_X+1, :] = cropped_I

    # same for segmentation
    seg_cell[_lb_Y-lb_Y+1 : -lb_Y+_ub_Y+1,
             _lb_X-lb_X+1 : -lb_X+_ub_X+1] = cropped_seg
    # keep only segmentation for cell i
    seg_cell[seg_cell .!= i] .= 0

    return I_cell, seg_cell
end