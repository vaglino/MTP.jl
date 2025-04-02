#= 
Background.jl contains methods to correct images illumination and subtract background from images
Stefano Travaglino, 2023, Zhu Lab, Georgia Tech
=#

module IlluminationCorrection

using TiffImages, Images, Statistics #, ImageFiltering, ImageMorphology
using PyCall
flt = pyimport("skimage.filters")
morph = pyimport("skimage.morphology")
cv2 = pyimport("cv2")
np = pyimport("numpy")

export average_illumination, illumination_correction, remove_outliers, illumination_correction_MTP, rescale_illumination_to_I


"""
Calculate average illumination of image stack
# Arguments
- `tifs_stack`: stack of tif image names to be used for illumination calculation
- `directory`: directory of tif images
- `channel`: channel to be used for illumination calculation (default = 1 for MTP, 2 for RICM)
# Returns
- `illumination`: average illumination of image stack 
"""
function average_illumination(tifs_stack; channel=1, directory="")
    MTP_all = [] 

    # first remove outliers for all images 
    for tif in tifs_stack
        img = Float64.(TiffImages.load(directory*tif))
        @show tif
        I = img[:,:,channel] .* 2^16 
        I_no_outl = remove_outliers(I)
        push!(MTP_all,I_no_outl)
    end

    # then take mean of all images 
    illumination = mean(MTP_all)
    # illumination = skexposure.rescale_intensity(illumination) .+ 1
    # finally, normalize illumination 0-1
    return illumination
end


"""
Correct illumination of image based on average illumination
# Arguments
- `I`: image to be corrected
- `illumination`: average illumination of image stack
# Returns
- `I_corrected`: corrected image
"""
function illumination_correction(I,illumination; darkI=zeros(size(I)))
      
    I_corrected = (I .- darkI) ./ illumination 

    # try later:
        # final_image = (corrected_image - min_intensity) / (illumination_profile - min_intensity)
    I_ = I .- darkI
    I_corrected = (I_ .- minimum(I_)) ./ (illumination .- minimum(I_))

    # subtract background value
    # background_value = mean(I_corrected) #+ std(I_corrected)
    # I_subtracted = I_corrected .- background_value

    # # remove background
    # I_subtracted = I_subtracted .- I_subtracted .* (I_subtracted .< 0)
end

function rescale_illumination_to_I(I,Illum, Illum_smooth)
    # Calculate the min and max for both images
    min_sample, max_sample = extrema(I)
    min_ref, max_ref = extrema(Illum_smooth)

    # Linearly scale the intensities of the background image
    rescaled_illum = min_sample .+ ((Illum .- min_ref) .* (max_sample - min_sample) ./ (max_ref - min_ref))
end

function illumination_correction_MTP(I,illumination,illumination_smooth)
    # remove image outliers and smooth image
    I_ = remove_outliers(I)
    # imshow(I_)
    I_smooth = imfilter(I_, Kernel.gaussian(50))
    
    filter_size = (31, 31)
    I_smooth = mapwindow(minimum, I_smooth, filter_size)
    # imshow(I_smooth)

    # rescale illumination to match image
    illum_resc = rescale_illumination_to_I(I_smooth,illumination, illumination_smooth)

    # correct image based on rescaled illumination
    I_corrected = illumination_correction(I_,illum_resc)
end

"""
Remove outliers from image based on neighborhood median
# Arguments
- `I`: image to be corrected
- `radius`: radius of disk kernel for median filter
- `threshold`: threshold for outlier detection
# Returns
- `_I`: image with outliers removed
"""
# function remove_outliers(I;radius=15,threshold=2) # :::TO DO::: make radius and threshold adaptive and add input argument for threshold values
#     # calculate median filtered image with disk kernel of given radius
#     median_filtered = flt.median(I,morph.disk(radius))
#     # imshow(median_filtered)

#     # outliers = I .> (median_filtered .+ threshold)
#     outliers_high = I .> (median_filtered .* (threshold))
#     outliers_low = I .< (median_filtered .* (1-0.3))
#     # imshow(outliers_high)
#     # imshow(outliers_low)

#     # combine outliers
#     outliers = Bool.(outliers_high .+ outliers_low)
#     # imshow(outliers)
#     # replace outliers with median filtered image
#     _I = deepcopy(I)
#     _I[outliers] = median_filtered[outliers]
#     return _I
# end
function remove_outliers(I::Matrix{Float64}; radius::Int=15, threshold::Float64=2.0)
    # Find the actual range of the input data
    min_val, max_val = extrema(I)
    
    # Scale to 8-bit range
    I_uint8 = UInt8.(round.((I .- min_val) ./ (max_val - min_val) .* 255))
    
    # Create a disk-shaped kernel
    kernel = morph.disk(radius)
    
    # Apply median filter using OpenCV
    median_filtered = cv2.medianBlur(I_uint8, 2*radius+1)
    
    # Convert back to Float64 and rescale to original range
    median_filtered_float = Float64.(median_filtered) ./ 255 .* (max_val - min_val) .+ min_val
    
    # Detect outliers (using broadcasting for efficiency)
    outliers = @. (I > median_filtered_float * threshold) | (I < median_filtered_float * (1 - 0.3))
    
    # Replace outliers
    _I = copy(I)
    @. _I[outliers] = median_filtered_float[outliers]
    
    return _I
end
end


##_______________________________________________________________________

module GaussianBackgroundSubtraction

using Distributions, Optim, Images
using LinearAlgebra

export fit_2d_gaussian, sample_gaussian, subtract_background, arr2tens


# Regularization term to ensure positive definiteness
const REGULARIZATION_TERM = 1e-6


# pack 2d gaussian variables
decode(x) = (μ = x[1:2], 

            σ = [x[3] x[4]; 
                 x[4] x[5]], 

            w =  x[6])

            # μ

function ensure_positive_definite(σ)
    # Add a small regularization term to the diagonal
    σ += REGULARIZATION_TERM * I
    return σ
end

function loss(I,p,I_ind)
    μ, σ, w = decode(p)
    σ = ensure_positive_definite(σ)
    try # σ might be infeasible so we have to handle this case

        tmp_dist = MvNormal(μ, σ)
        y = map( pxl -> w * pdf(tmp_dist,pxl) , I_ind) 
        MSE = sum( (I.-y).^2 )
    catch
        sum(I)
    end
end

function fit_2d_gaussian(I)
    I_ind =[ [i,j] for i in 1:size(I,1), j in 1:size(I,2)]
    _loss = (p) -> loss(I,p,I_ind)
 
    m = n = size(I,1)
    p_init = [m/2, m/2, m, 0, m, m]

    res = optimize(_loss, p_init , LBFGS() )#; optimizer_o = Optim.Options(show_trace = true)) # uncomment for verbose output
    μ, σ, w = decode(res.minimizer)
    res.minimizer
end

function sample_gaussian(p,I)
    μ, σ, w = decode(p) # unpack multigaussian parameters

    σ = ensure_positive_definite(σ)
    # Check for positive definiteness
    if det(σ) <= 0
        throw(ArgumentError("The covariance matrix is not positive definite"))
    end

    # matrix of indexes
    I_ind =[ [i,j] for i in 1:size(I,1), j in 1:size(I,2)]

    # 2d gaussian distribution
    tmp_dist = MvNormal(μ, σ) 
    #sample gaussian at indexes
    y = map( pxl -> w * pdf(tmp_dist,pxl) , I_ind)
    
    return y 
end

"""
Subtract background from image using 2d gaussian fit.
# Arguments
- `I`: image to be background subtracted
- `channels`: channels to be background subtracted
- `sf`: subsampling factor, to speed up fit
# Returns
- `I_subtracted`: background subtracted image
- `gauss`: background gaussian
"""
function subtract_background(I; channels = nothing, sf = 30)

    if channels === nothing
        channels = 1:size(I,3)
    end
    
    gauss = []
    I_subtracted = []

    # for k in 3:3
    for k in channels

        Iₖ = deepcopy(I[:,:,k])
        # imshow(Iₖ)
        
        # subsample image to speed up fit
        Iₖ_small = Iₖ[1:sf:end,1:sf:end]

        # fit 2d gaussian
        pₖ = fit_2d_gaussian(Iₖ_small) 
        gaussₖ_small = sample_gaussian(pₖ,Iₖ_small)
        
        # upsample 2d gaussian to original size
        gaussₖ = imresize(gaussₖ_small, size(Iₖ))

        # subtract background gaussian from image
        Iₖ_subtracted = Iₖ - gaussₖ
        Iₖ_subtracted = Iₖ_subtracted .- minimum(Iₖ_subtracted)
        # imshow(Iₖ_subtracted)

        push!(I_subtracted,Iₖ_subtracted)
        push!(gauss,gaussₖ)
    end

    return I_subtracted, gauss
end

function arr2tens(arr)
    I = arr[1]
    for i in 2:length(arr)
        I = cat(I,arr[i],dims=3)
    end
    return I
end

end

##_______________________________________________________________________

module RollingBallSubtraction

using PyCall

skexposure = pyimport("skimage.exposure")
restoration = pyimport("skimage.restoration")

export rolling_ball

restoration = pyimport("skimage.restoration")

"""
Background subtraction using rolling ball method
# Arguments
- `I`: image to be background subtracted
- `radius`: radius of rolling ball
- `channels`: channels to be background subtracted
# Returns
- `backgrounds`: background subtracted images
"""
function rolling_ball(I;radius=200,channels=nothing) # :::TO DO::: change call for this function to include channels=1 argument
    if channels === nothing
        channels = 1:size(I,3)
    end

    backgrounds = []
    for i in channels
        _I = I[:,:,i]
        backgroundᵢ = restoration.rolling_ball(_I,radius=radius)
        push!(backgrounds, backgroundᵢ)
    end
    return backgrounds
end

end


using Images, ImageFiltering, ImageMorphology


function create_disk(radius)
    r2 = radius * radius
    disk_element = [sqrt(I[1]^2 + I[2]^2) <= radius for I in CartesianIndices((-radius:radius, -radius:radius))]
    return disk_element
end

function remove_outliers2(I; radius=30, threshold=2)
    # Create a disk-shaped structuring element
    disk_element = create_disk(radius)

    # Calculate median-filtered image with disk kernel of the given radius
    median_filtered =  mapwindow(median!, I, (15,15))
    
    # Identify high and low outliers
    outliers_high = I .> (median_filtered .* threshold)
    outliers_low = I .< (median_filtered .* (1 - 0.3))
    
    # Combine outliers
    outliers = outliers_high .| outliers_low
    
    # Replace outliers with median-filtered image
    I_filtered = deepcopy(I)
    I_filtered[outliers] = median_filtered[outliers]
    return I_filtered
end