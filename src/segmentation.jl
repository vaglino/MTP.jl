# this file contains all segmentation functions for every channel
# add as needed, maintaining the same conventions

using PyCall
morph = pyimport("skimage.morphology")
feat = pyimport("skimage.feature")
flt = pyimport("skimage.filters")
ndimage = pyimport("scipy.ndimage")
skseg = pyimport("skimage.segmentation")
skexposure = pyimport("skimage.exposure")


function segment(I; layer="DIC")

    if layer == "RICM"
        segment_RICM(I)
    elseif layer == "DIC"
        segment_DIC(I)
    elseif layer == "TIRF"
        segment_TIRF(I)
    else
        @show "unknown image type, will segment assuming TIRF"
        segment_TIRF(I)
    end

end

function segment_DIC(I)
    plotting = false
    #equalize contrast
    I_eq = adjust_histogram(I, AdaptiveEqualization(nbins = 256, rblocks = 8, cblocks = 8, clip = 0.1))
    plotting && imshow(I_eq)

    #gaussian smooth to remove noise
    I_eq_smooth = flt.gaussian(I_eq,sigma=3)
    plotting && imshow(I_eq_smooth)
    
# #____________________
#     #gaussian smooth to remove noise
#     I_smooth = flt.gaussian(I,sigma=5)
#     imshow(I_smooth)

#     #equalize contrast
#     I_eq_smooth = adjust_histogram(I_smooth, AdaptiveEqualization(nbins = 256, rblocks = 8, cblocks = 8, clip = 0.))
#     imshow(I_eq_smooth)
# #______________
    #edge detection
    # edges = flt.sobel(I_eq_smooth)
    edges = feat.canny(I_eq_smooth, sigma=3)
    plotting && imshow(edges)

    #img dilation
    dilated = morph.binary_dilation(edges,morph.disk(5))
    plotting && imshow(dilated)

    #fill holes
    filled =ndimage.binary_fill_holes(dilated)#,structure=morph.diamond(5))
    plotting && imshow(filled)

    #erode twice
    eroded = ndimage.binary_erosion(filled,morph.disk(5) )
    eroded2 = ndimage.binary_erosion(eroded,morph.disk(5) )

    plotting && imshow(eroded2)
    bw = morph.remove_small_objects(eroded2,500)
    plotting && imshow(bw)

    bw2 = skseg.clear_border(bw)

    #watershed
    segments = watershed_segmentation(bw2)
    # readline()
    return segments, bw2
end


function segment_RICM(I; seeds = nothing)
    fuzz_factor  = 0.9
    I_inverted = 1 .- I
    thres = otsu_threshold(I_inverted)
    # bw = (I_inverted) .> thres
    bw = (I_inverted) .>  thres .* fuzz_factor

    Iopen = opening(bw)
    Iopen = opening(Iopen)
    Iclose = closing(Iopen)

    Ifilled = ndimage.binary_fill_holes(Iclose)
 
    bw3 = morph.remove_small_objects(Ifilled,2000)
    # imshow(bw3)
    bw4 = skseg.clear_border(bw3)
    bw4 = dilate(dilate(bw4))
    
    segments = watershed_segmentation(bw4)
    # segments_labels = labels_map(segments) .* bw4
    # segments_no_border = skseg.clear_border(segments_labels)
    
    return segments, bw4
    # return segments_no_border, bw4
end



function segment_RICM_single_cell(I; seeds = nothing)
    plotting = false
    #-------------------
    I_rsc = skexposure.rescale_intensity(I)#,in_range=(0.7*minimum(I), 1*maximum(I)))
    # I_eq = adjust_histogram(I, AdaptiveEqualization(nbins = 256, rblocks = 8, cblocks = 8, clip = 0.1))
    # imshow(I_rsc)
    I_inverted = 1 .- I_rsc
    # imshow(I_inverted)
    #------------------------------
    # I_eq = skexposure.equalize_hist(I_rsc)
    # I_inverted = 1 .- I_eq

    #gaussian smooth to remove noise
    I_smooth = flt.gaussian(I_inverted,sigma=1) #1
    plotting && imshow(I_smooth)
    I_inverted = I_smooth

    # I_inverted = skexposure.equalize_hist(I_inverted)
    # plotting &&  imshow(I_inverted)

    thres = otsu_threshold(I_inverted)
    bw = (I_inverted) .> thres
    # plotting && imshow(bw)

    Iopen = opening(bw)
    Iopen = opening(Iopen)
    #imshow(Iopen)
    Iclose = closing(Iopen)
    # imshow(Iclose)
    Ifilled = ndimage.binary_fill_holes(Iclose)
    # Ifilled = dilate(dilate(Ifilled))
    plotting &&  imshow(Ifilled)

    bw3 = morph.remove_small_objects(Ifilled,3/px_area) #3
    plotting &&  imshow(bw3)


    bw5 = morph.remove_small_objects(bw3,3/px_area)
    bw6 = morph.binary_dilation(bw5,morph.disk(10)) #10

    plotting &&  imshow(bw6)

    segments = watershed_segmentation(bw6)
    # @show segments
    

    labels = labels_map(segments)
    labels = labels .* bw6
    if plotting
        overlay = plot_segmentation_overlay(I_rsc,labels)
        imshow(overlay)
    end
    # plotting && imshow(labels)
    # if maximum(labels) > 1
    #     segments_clear = skseg.clear_border(labels)
    # else
    #     segments_clear = labels
    # end
    segments_clear = skseg.clear_border(labels)
    plotting && imshow(segments_clear)
    ricm_segmented = segments_clear .!= 0
    all_segments = labels .!=0
    # ricm_segmented = morph.binary_dilation(ricm_segmented,morph.disk(5)) #comment
    plotting && imshow(ricm_segmented)

    return segments, ricm_segmented, all_segments
    # return segments_no_border, bw4
end


function watershed_segmentation(bw)
    plotting = false

    dist =  - distance_transform(feature_transform(.!bw))
    dist = flt.gaussian(dist,sigma=5)
    plotting && imshow(dist)

  
    # mins = feat.peak_local_max(-dist,min_distance=50,threshold_rel=0.1)
    mins = feat.peak_local_max(-dist,min_distance=50,threshold_rel=0.1,exclude_border=false)
    mins = mins .+ 1 #sccount for python 0 indexing vs. julia 1 indexing
    # @show mins
    mins = [CartesianIndex(mins[i,1],mins[i,2]) for i in 1:size(mins,1)]
    # @show mins

    markers = zeros(size(dist))
    markers[mins] .= 1

    #include highest nonzero marker 
    edges = [dist[1,:],]


    markers = morph.binary_dilation(markers,morph.diamond(5))
    plotting &&  imshow(markers)
   
    # markers_dist = label_components(dist .< -50.)
    markers = Bool.(markers)
    # markers = morph.binary_dilation(markers,morph.diamond(5))
    # imshow(markers)
    markers = label_components(markers)
    
    segments = watershed(dist, markers) 
end



using IndirectArrays
function plot_segmented(I_seg, I_bw)
    labels = labels_map(I_seg)
    # labels = I_seg
    colored_labels = IndirectArray(labels, distinguishable_colors(maximum(labels)))
    masked_colored_labels = colored_labels .* (I_bw)
    # mosaic(I_corrected, colored_labels, masked_colored_labels; nrow=1)
    return masked_colored_labels
end

skcolor = pyimport("skimage.color")
function plot_segmentation_overlay(I,labels)
    I_rsc = skexposure.rescale_intensity(I) 
    I_overlay = Float64.(skcolor.label2rgb(labels,I_rsc,bg_label=0))
    I_over_perm = permutedims(I_overlay, (3, 1, 2))
    overlay = colorview(RGB,I_over_perm)
    # imshow(overlay)
    return overlay
end