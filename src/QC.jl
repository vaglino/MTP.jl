
# function QC(cells; condition_name="Condition X")

#     correct = []
#     for (i, cell) in enumerate(cells)
#         overlay_MTP = mosaic_cropped_cells([cell];channel=1)
#         # set NAN values to 0 in overlay_MTP
#         @show sum(isnan.(overlay_MTP[1]))
#         overlay_MTP[1][isnan.(overlay_MTP[1])] .= 0
        

#         plt1 = plot(overlay_MTP[1], title="MTP")
#         overlay_RICM = mosaic_cropped_cells([cell];channel=2)
#         plt2 = plot(overlay_RICM[1], title="RICM")
#         overlay_DIC = mosaic_cropped_cells([cell];channel=3)
#         plt3 = plot(overlay_DIC[1], title="DIC")
#         plt4 = plot(Gray.( skexposure.rescale_intensity( cell.I_fixed[:,:,3] ) ), title="DIC (fixed crop)")

#         plot_size = (1000, 1000)
#         overall_title = string(condition_name, " -- cell # ", i, " / ", length(cells))
#         plt = plot(plt1, plt2, plt3, plt4, title=["MTP" "RICM" "DIC" "DIC (fixed crop)"], size=plot_size,plot_title=overall_title)
#         display(plt)

#         println("Correct (0 or 1)")

#         input = nothing
#         while input != 0 && input != 1
#             try
#                 input_str = readline()
#                 input = parse(Int64, input_str)
#                 if input != 0 && input != 1
#                     println("Invalid input. Please enter either 0 or 1:")
#                 end
#             catch e
#                 if isa(e, InterruptException)
#                     println("Exiting QC function.")
#                     return correct
#                 else
#                     println("Invalid input. Please enter either 0 or 1:")
#                 end
#             end
#         end

#         push!(correct,input)
#     end
#     return correct
# end
    
function QC(cells; condition_name="Condition X")
    correct = []
    
    # Pre-process all channels for all cells at once
    all_channels = Dict(i => Dict() for i in 1:length(cells))
    for (i, cell) in enumerate(cells)
        for channel in 1:3
            all_channels[i][channel] = mosaic_cropped_cells([cell]; channel=channel)[1]
        end
    end
    
    for (i, cell) in enumerate(cells)
        # Use pre-processed images
        plt1 = plot(all_channels[i][1], title="MTP")
        plt2 = plot(all_channels[i][2], title="RICM")
        plt3 = plot(all_channels[i][3], title="DIC")
        plt4 = plot(Gray.(skexposure.rescale_intensity(cell.I_fixed[:,:,3])), title="DIC (fixed crop)")
        
        plot_size = (1000, 1000)
        overall_title = string(condition_name, " -- cell # ", i, " / ", length(cells))
        plt = plot(plt1, plt2, plt3, plt4, title=["MTP" "RICM" "DIC" "DIC (fixed crop)"], size=plot_size, plot_title=overall_title)
        display(plt)
        
        println("Correct (0 or 1)")

        input = nothing
        while input != 0 && input != 1
            try
                input_str = readline()
                input = parse(Int64, input_str)
                if input != 0 && input != 1
                    println("Invalid input. Please enter either 0 or 1:")
                end
            catch e
                if isa(e, InterruptException)
                    println("Exiting QC function.")
                    return correct
                else
                    println("Invalid input. Please enter either 0 or 1:")
                end
            end
        end

        push!(correct,input)
    end
    
    return correct
end
# function QC(cells)

#     correct = []
#     for (i, cell) in enumerate(cells)
#         overlay_MTP = mosaic_cropped_cells([cell];channel=1)
#         plt1 = plot(overlay_MTP[1])
#         overlay_RICM = mosaic_cropped_cells([cell];channel=2)
#         plt2 = plot(overlay_RICM[1])
#         overlay_DIC = mosaic_cropped_cells([cell];channel=3)
#         plt3 = plot(overlay_DIC[1])
#         plt4 = plot(Gray.( skexposure.rescale_intensity( cell.I_fixed[:,:,3] ) ) )

#         plot_size = (1000, 1000)
#         plt = plot(plt1, plt2, plt3, plt4, title=["MTP" "RICM" "DIC" "DIC (fixed crop)"], size=plot_size)
#         display(plt)
#         # plt = plot(cell)
#         # display(plt)
#         println("cell # ", i, " Correct (0 or 1)")
#         println("Correct (0 or 1)")
#         # input = readline()
#         # try 
#         #     input = parse(Int64, input) 
#         # catch
#         #     println("re enter input")
#         #     input = readline()
#         #     input = parse(Int64, input) 
#         # end
#         input = nothing
#         while input != 0 && input != 1
#             try
#                 input_str = readline()
#                 input = parse(Int64, input_str)
#                 if input != 0 && input != 1
#                     println("Invalid input. Please enter either 0 or 1:")
#                 end
#             catch
#                 println("Invalid input. Please enter either 0 or 1:")
#             end
#         end
#         @show input

#         # if input == 0 | input == 1
#         #     correct = input
#         # else
#         #     println("re enter input")
#         #     input = readline()
#         #     correct = input
#         # end
#         push!(correct,input)
#         # close(img)
#     end
#     return correct
# end

# # QC_res = QC(cells)
# # QC_res[QC_res .==2] .= 0
# # QC_res
