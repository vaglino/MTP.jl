# Molecular Tension Probe (MTP) Analysis

This Julia library provides a complete workflow for analyzing Molecular Tension Probe microscopy data to quantify mechanical forces exerted by cells. The library processes multi-channel microscopy images (MTP, RICM, and DIC), segments individual cells, and extracts quantitative tension measurements.

## Overview

The MTP Analysis library analyzes images from molecular tension probe experiments, which use fluorescent sensors to measure cellular traction forces. It corrects for uneven illumination, segments cells, extracts quantitative measurements, and provides visualization tools for the results.

## Basic Usage

```julia
# Import necessary modules
using Revise
using Images, ImageIO
using Statistics, Plots
using TiffImages, JLD2

# Add path to source files
src = "/path/to/source/files/"
includet(src*"segmentation.jl")
includet(src*"background.jl")
includet(src*"image_processing.jl")
includet(src*"QC.jl")
includet(src*"analysis_functions.jl")
includet(src*"scripting_functions.jl")

using .IlluminationCorrection

# Set pixel size parameter
px_size = 0.1083333 # μm
px_area = px_size^2 # μm²

# Define data and result directories
data_dir = "/path/to/data/directory/"
res_dir = "/path/to/results/directory/"

# Get conditions from subdirectories
conds = get_conditions(data_dir)
print(keys(conds))

# Calculate background illumination
illumination_MTP = average_illumination(conds["back"][2]; channel=1, directory=conds["back"][1])
imshow(illumination_MTP)
illumination_RICM = average_illumination(conds["back"][2]; channel=2, directory=conds["back"][1])
imshow(illumination_RICM)

# save average background illumination files
save_object(res_dir*"illumination_MTP.jld2", illumination_MTP)
save_object(res_dir*"illumination_RICM.jld2", illumination_RICM)

# Process data
results = process_data(data_dir; illumination = illumination_MTP)

# Perform quality control
gr()
QC_results, ΣFIs_all, MFIs_all, SAs_all = perform_QC(results)

# Plot results
plotly()
plot_results(QC_results)

# Save results
save_conditions_to_csv(QC_results, res_dir*"experiment_name")
save_object(res_dir*"experiment_name_cells.jld2", results)

#re-load objects if necessary, without having to repeat the analysis from scratch
# results = load_object(res_dir*"experiment_name_cells.jld2")
# illumination_MTP = load_object(res_dir*"illumination_MTP.jld2")
# illumination_RICM = load_object(res_dir*"illumination_RICM.jld2")
```

## Installation

### Prerequisites
- Julia 1.6+
- Required Julia packages (listed in `library list.txt`)
- Python with scikit-image, scipy, numpy, and opencv

### Required Packages
```julia
using Pkg
Pkg.add([
    "CSV", "DataFrames", "Distributions", "ImageBinarization", 
    "ImageFiltering", "ImageIO", "ImageMorphology", "ImageSegmentation", 
    "ImageView", "Images", "IndirectArrays", "JLD2", "LinearAlgebra", 
    "Optim", "Plots", "PyCall", "Random", "Revise", "StaticArrays", 
    "Statistics", "StatsPlots", "TiffImages"
])
```

### Python Dependencies
```bash
pip install scikit-image scipy numpy opencv-python
```

## Data Directory Structure
The library expects your data to be organized as follows:
```
data_dir/
    |-> Condition1/
        |-> image_1.tif
        ...
    |-> Condition2/
        |-> image_1.tif
        ...
    |-> back/  # Background reference images
        |-> back_1.tif
        ...
```

## Analysis Pipeline

### 1. Illumination Correction
The library first corrects for uneven illumination using background reference images:

```julia
# Find MTP and RICM background illumination patterns
illumination_MTP = average_illumination(conds["back"][2]; channel=1, directory=conds["back"][1])
illumination_RICM = average_illumination(conds["back"][2]; channel=2, directory=conds["back"][1])
```

This step creates average illumination profiles from background images to normalize all experimental images.

### 2. DIC Image Segmentation
Cells are segmented from the DIC channel using a combination of filters, edge detection, and watershed segmentation. This process identifies individual cells in the field of view.

### 3. Cell Object Creation
The library creates a `Cell` struct for each segmented cell, containing:
- Cropped multi-channel images (`I`)
- Segmentation masks (`seg`)
- Fixed-size representations for display (`I_fixed`, `seg_fixed`)
- Coordinate information

This Cell object is key for advanced usage, as it allows to peer into each single cell separately. For example, after completing the analysis, one may want to check the images and segmentations of a specific cell.

```julia
# get the 10th cell of a specific condition
cell_10 = results["condition_name"][10]

# plot MTP channel of the 10th cell
imshow(cell_10.I[:,:,1])
# plot MTP channel of the 10th cell
imshow(cell_10.I[:,:,2])
# plot MTP channel of the 10th cell
imshow(cell_10.I[:,:,3])

# plot MTP channel of the 10th cell with fixed crop ratio (useful for representative images selection)
imshow(cell_10.I_fixed[:,:,1])
```

### 4. RICM Segmentation
For each cell, the RICM channel is segmented to identify the cell-substrate contact area. This provides a mask for measuring tension metrics only in areas where the cell is in contact with the substrate.

### 5. Manual Quality Control
The library provides an interactive interface for reviewing cell segmentations and accepting/rejecting each cell:

```julia
# Perform manual QC on all conditions
QC_results, ΣFIs_all, MFIs_all, SAs_all = perform_QC(results)
```

This step ensures that only properly segmented cells are included in the final analysis.

### 6. Analysis and Visualization
After QC, the library calculates metrics and visualizes the results:
- Sum of Fluorescence Intensity (ΣFI)
- Mean Fluorescence Intensity (MFI)
- Spread Area (SA)

## Advanced Usage

### Customizing Segmentation Parameters

You can modify the segmentation parameters to improve cell detection for your specific images:

```julia
# Modify DIC segmentation in segmentation.jl
function segment_DIC(I)
    # Customize contrast adjustment
    I_eq = adjust_histogram(I, AdaptiveEqualization(nbins = 256, rblocks = 10, cblocks = 10, clip = 0.15))
    
    # Customize Gaussian smoothing
    I_eq_smooth = flt.gaussian(I_eq, sigma=4)  # Increase from default 3 for more smoothing
    
    # Customize edge detection
    edges = feat.canny(I_eq_smooth, sigma=5)  # Increase from default 4 for more robust edges
    
    # Customize dilation
    dilated = morph.binary_dilation(edges, morph.disk(7))  # Increase from default 5 for larger dilation
    
    # Customize small object removal
    bw = morph.remove_small_objects(eroded, 3000)  # Increase from default 2000 for larger minimum cell size
    
    # Rest of the function remains the same...
end
```

### Changing the Crop Ratio Around Single Cells

You can modify the crop ratio to adjust how much area around a cell is included in the `Cell` object:

```julia
# In image_processing.jl, modify the crop_cell function
function crop_cell(I, seg, i; crop_ratio=2.0)  # Increase from default 1.5 for larger crop area
    # Implementation remains the same...
end
```

Alternatively, directly modify the call to crop_cell:

```julia
# When calling crop_cell
I_cell, segmented = crop_cell(I, seg, i, crop_ratio=2.5)  # Custom crop ratio
```

### Changing Fixed Size Crop

You can modify the fixed size for standardized cell visualization:

```julia
# In image_processing.jl, modify the crop_cell_fixed_size function
function crop_cell_fixed_size(I, seg, i)
    # Change the crop size (in micrometers)
    crop_size = round(Int, 20 / px_size)  # Change from default 17.5 μm
    
    # Rest of the function remains the same...
end
```

### Customizing RICM Segmentation

The RICM channel segmentation identifies the cell-substrate contact area and is crucial for accurate tension measurements. You can customize the RICM segmentation parameters in the `segment_RICM_single_cell` function:

```julia
# In segmentation.jl, modify the segment_RICM_single_cell function
function segment_RICM_single_cell(I; seeds = nothing)
    # Customize Gaussian smoothing (increase for more noise reduction)
    I_smooth = flt.gaussian(I_inverted, sigma=2)  # Change from default sigma=1
    
    # Use a different thresholding algorithm
    # algo = Otsu()  # Default
    algo = Yen()     # Alternative thresholding algorithm
    bw = binarize(I_inverted, algo)
    
    # Adjust morphological operations
    # Apply more opening operations for smoother boundaries
    Iopen = opening(opening(opening(bw)))
    
    # Modify small object removal size threshold
    bw3 = morph.remove_small_objects(Ifilled, 2/px_area)  # Increase from default 1/px_area
    
    # Adjust dilation radius (increase for larger contact area detection)
    bw6 = morph.binary_dilation(bw5, morph.disk(15))  # Increase from default 10
    
    # Rest of the function remains the same...
end
```

Fine-tuning these parameters can help improve RICM segmentation for:
- Cells with different morphologies
- Images with varying contrast or noise levels
- Different substrate compositions affecting RICM contrast
- Experiments with different cell types that have unique adhesion patterns

## Source Code Overview

### segmentation.jl
Contains cell segmentation functions:
- `segment_DIC`: Segments cells in DIC images
- `segment_RICM_single_cell`: Segments contact area in RICM images
- `watershed_segmentation`: Separates touching cells

### background.jl
Provides illumination correction:
- `average_illumination`: Creates average illumination profile
- `illumination_correction`: Corrects images using the profile
- `remove_outliers`: Removes noise and artifacts

### image_processing.jl
Processes images and creates cell objects:
- `Cell` structure definition
- `files2cells`: Processes multiple image files
- `img2cells`: Extracts cells from a single image
- `add_RICM_segmentation_to_cell`: Adds contact area segmentation

### analysis_functions.jl
Quantifies cellular tension:
- `MFI`: Calculates mean fluorescence intensity
- `ΣFI`: Calculates sum of fluorescence intensity
- `spread_area`: Calculates cell spread area

### QC.jl
Provides manual quality control:
- `QC`: Interactive function for reviewing segmentations

### scripting_functions.jl
High-level workflow functions:
- `get_conditions`: Identifies experimental conditions
- `process_data`: Processes all conditions
- `perform_QC`: Runs quality control
- `plot_results`: Visualizes analysis
- `save_conditions_to_csv`: Exports results

