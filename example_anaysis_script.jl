using Revise
using Images, ImageView, ImageIO
using ImageSegmentation
using Statistics, Random, Plots
using PyCall
using TiffImages
using JLD2
using StatsPlots

src = "/Users/stefano/Documents/ZhuLab/MTP_v2/src/"
includet(src*"segmentation.jl")
includet(src*"background.jl")
includet(src*"image_processing.jl")
includet(src*"QC.jl")
includet(src*"analysis_functions.jl")
includet(src*"scripting_functions.jl")

using .IlluminationCorrection, .GaussianBackgroundSubtraction#, .RollingBallSubtraction

# SET PARAMETERS
px_size = 0.1083333 #um
px_area = px_size^2 #um²

# DATA DIRECTORY (set up data directory as follows)
# data_dir
#       |-> Condition 1
#          |-> image_1.tif
#          ...
#          ...
#          |-> image_N.tif
#       ...
#       ...
#       |-> Condition N


#=_______________________________________________________________________________________________________________________
# IgG ITTYF BSA-NIP gradient 4.7 pN (20250321)
_______________________________________________________________________________________________________________________=#

res_dir = "/Volumes/labs4/zhu-lab/Stefano/MTP/20250321_IgG ITT2F BSA-NIP/20250321_IgG ITT2F BSA-NIP/analyzed_cells_4.7pN/"

data_dir = "/Volumes/labs4/zhu-lab/Stefano/MTP/20250321_IgG ITT2F BSA-NIP/20250321_IgG ITT2F BSA-NIP/4.7pN/"
# data_dir = "Y:/zhu-lab/Amir/MTP/20250114/20250114/19 pN/"


# check conditions for each subdirectory
conds = get_conditions(data_dir)
print(keys(conds))

# # Find MTP background illumination by averaging the BSA surfaces
illumination_MTP = average_illumination(conds["back"][2];channel=1, directory=conds["back"][1])
imshow(illumination_MTP)
# Find RICM background illumination by averaging the BSA surfaces
illumination_RICM = average_illumination(conds["back"][2];channel=2, directory=conds["back"][1])
imshow(illumination_RICM)


@time results = process_data(data_dir; illumination = illumination_MTP)
# imshow(results["0.01"][47].I[:,:,1])
gr()
QC_results, ΣFIs_all, MFIs_all, SAs_all = perform_QC(results)

plotly()
# plot_results(ΣFIs_all, MFIs_all, SAs_all)
plot_results(QC_results)

#-----------------------------------------------------------------------
## SAVE ANALYSIS
save_conditions_to_csv(QC_results, res_dir*"20250321_IgG ITT2F BSA-NIP_4.7pN")
save_object(res_dir*"20250321_IgG ITT2F BSA-NIP_4.7pN_cells.jld2", results)
save_object(res_dir*"illumination_MTP_4.7pN.jld2", illumination_MTP)
save_object(res_dir*"illumination_RICM_4.7pN.jld2", illumination_RICM)
results = load_object(res_dir*"20250321_IgG ITT2F BSA-NIP_4.7pN_cells.jld2")
illumination_MTP = load_object(res_dir*"illumination_MTP_4.7pN.jld2")
illumination_RICM = load_object(res_dir*"illumination_RICM_4.7pN.jld2")
#-----------------------------------------------------------------------


#=_______________________________________________________________________________________________________________________
# IgG ITTYF BSA-NIP gradient 12 pN (20250321)
_______________________________________________________________________________________________________________________=#
res_dir = "/Volumes/labs4/zhu-lab/Stefano/MTP/20250321_IgG ITT2F BSA-NIP/20250321_IgG ITT2F BSA-NIP/analyzed_cells_12pN/"
# data_dir = "Y:/zhu-lab/Amir/MTP/20250114/20250114/analyzed_cells/"


data_dir = "/Volumes/labs4/zhu-lab/Stefano/MTP/20250321_IgG ITT2F BSA-NIP/20250321_IgG ITT2F BSA-NIP/12pN/"
# data_dir = "Y:/zhu-lab/Amir/MTP/20250114/20250114/19 pN/"


# check conditions for each subdirectory
conds = get_conditions(data_dir)
print(keys(conds))

# # # Find MTP background illumination by averaging the BSA surfaces
illumination_MTP = average_illumination(conds["back"][2];channel=1, directory=conds["back"][1])
imshow(illumination_MTP)
# Find RICM background illumination by averaging the BSA surfaces
illumination_RICM = average_illumination(conds["back"][2];channel=2, directory=conds["back"][1])
imshow(illumination_RICM)


@time results_12pN = process_data(data_dir; illumination = illumination_MTP)
# imshow(results["0.01"][47].I[:,:,1])
gr()
QC_results, ΣFIs_all, MFIs_all, SAs_all = perform_QC(results_12pN)

plotly()
# plot_results(ΣFIs_all, MFIs_all, SAs_all)
plot_results(QC_results)

#-----------------------------------------------------------------------
## SAVE ANALYSIS
save_conditions_to_csv(QC_results, res_dir*"20250321_IgG ITT2F BSA-NIP_12pN")
save_object(res_dir*"20250321_IgG ITT2F BSA-NIP_12pN_cells.jld2", results_12pN)
save_object(res_dir*"illumination_MTP_12pN.jld2", illumination_MTP)
save_object(res_dir*"illumination_RICM_12pN.jld2", illumination_RICM)
results_12pN = load_object(res_dir*"20250321_IgG ITT2F BSA-NIP_12pN_cells.jld2")

#-----------------------------------------------------------------------


#=_______________________________________________________________________________________________________________________
# IgG ITTYF BSA-NIP gradient 19 pN (20250321)
_______________________________________________________________________________________________________________________=#
res_dir = "/Volumes/labs4/zhu-lab/Stefano/MTP/20250321_IgG ITT2F BSA-NIP/20250321_IgG ITT2F BSA-NIP/analyzed_cells_19pN/"
# data_dir = "Y:/zhu-lab/Amir/MTP/20250114/20250114/analyzed_cells/"


data_dir = "/Volumes/labs4/zhu-lab/Stefano/MTP/20250321_IgG ITT2F BSA-NIP/20250321_IgG ITT2F BSA-NIP/19pN/"
# data_dir = "Y:/zhu-lab/Amir/MTP/20250114/20250114/19 pN/"


# check conditions for each subdirectory
conds = get_conditions(data_dir)
print(keys(conds))

# # Find MTP background illumination by averaging the BSA surfaces
illumination_MTP = average_illumination(conds["back"][2];channel=1, directory=conds["back"][1])
imshow(illumination_MTP)
# Find RICM background illumination by averaging the BSA surfaces
illumination_RICM = average_illumination(conds["back"][2];channel=2, directory=conds["back"][1])
imshow(illumination_RICM)


@time results = process_data(data_dir; illumination = illumination_MTP)
# imshow(results["0.01"][47].I[:,:,1])
gr()
QC_results, ΣFIs_all, MFIs_all, SAs_all = perform_QC(results)

plotly()
# plot_results(ΣFIs_all, MFIs_all, SAs_all)
plot_results(QC_results)

#-----------------------------------------------------------------------
## SAVE ANALYSIS
save_conditions_to_csv(QC_results, res_dir*"20250321_IgG ITT2F BSA-NIP_19pN")
save_object(res_dir*"20250321_IgG ITT2F BSA-NIP_19pN_cells.jld2", results)
save_object(res_dir*"illumination_MTP_19pN.jld2", illumination_MTP)
save_object(res_dir*"illumination_RICM_19pN.jld2", illumination_RICM)
illumination_MTP = load_object(res_dir*"illumination_MTP_19pN.jld2")
illumination_RICM = load_object(res_dir*"illumination_RICM_19pN.jld2")

results_003 = Dict("0.03" => results["0.03"])
QC_results, ΣFIs_all, MFIs_all, SAs_all = perform_QC(results_003)
save_conditions_to_csv(QC_results, res_dir*"20250321_IgG ITT2F BSA-NIP_19pN_0_03_only")
results["0.03"] = results_003["0.03"]
results
plotly()
plot_results(QC_results)