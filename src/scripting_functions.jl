"""
Prints the conditions names found in the data directory
# Arguments
- `data_dir`: the path to the data directory
# Returns
- `conditions_dirs`: a dictionary with the condition names as keys and the corresponding directories and files as values
"""
function get_conditions(data_dir)
    subdirs = [d for d in readdir(data_dir) if isdir(joinpath(data_dir, d))]
    conditions_dirs = Dict()
    for subdir in subdirs
        dir = joinpath(data_dir, subdir *"/")
        files = readdir(dir)
        tifs = filter(file -> endswith(file, ".tif"), files)
        conditions_dirs[subdir] = (dir, tifs)
    end
    return conditions_dirs
end

"""
Processes a single condition by running files2cells and flattening the resulting list of cell objects
# Arguments
- `dir`: the path to the condition directory
- `files`: the files in the condition directory
- `illumination`: average illumination of image stack
# Returns
- `cells_flat`: a list of all cell objects in the condition
"""
function process_condition(dir, files, illumination)
    cells = files2cells(dir, files; illumination=illumination)
    ImageView.closeall()
    cells_flat = reduce(vcat, cells)
    return cells_flat
end

function process_QC(cells; condition_name="Condition X")
    QC_res = QC(cells)
    ΣFIs, MFIs, SAs = summarize_cells(cells, QC_res)
    return QC_res, ΣFIs, MFIs, SAs
end

"""
Processes all conditions in the data directory
# Arguments
- `data_dir`: the path to the data directory
- `illumination`: average illumination of image stack
# Returns
- `results`: a dictionary with the condition names as keys and the corresponding cell objects as values
"""
function process_data(data_dir; illumination=illumination_MTP)

    # Get conditions from subdirectories
    conditions_dirs = get_conditions(data_dir)

    # Process all conditions
    results = Dict()
    for (name, (dir, files)) in conditions_dirs
        cells = process_condition(dir, files, illumination)
        results[name] = cells
    end
    return results
end

   # Plotting code
using StatsPlots
"""
Performing manua QC for the results of the analysis of all conditions in the data directory
# Arguments
- `results`: a dictionary with the condition names as keys and the corresponding cell objects as values
# Returns
- `QC_results`: a dictionary with the condition names as keys and the corresponding QC results as values
- `ΣFIs_all`: a vector with the ΣFIs of all cells
- `MFIs_all`: a vector with the MFIs of all cells
- `SAs_all`: a vector with the SAs of all cells
"""
function perform_QC(results)
    # Perform QC for all conditions and extract the processed values
    QC_results = Dict()
    for name in keys(results)
        cells = results[name]
        QC_res = QC(cells, condition_name = name) # Call QC with the cells and the condition name
        ΣFIs, MFIs, SAs = summarize_cells(cells, QC_res)
        QC_results[name] = (QC_res, ΣFIs, MFIs, SAs)
    end

    ΣFIs_all = [QC_results[name][2] for name in keys(QC_results)]
    MFIs_all = [QC_results[name][3] for name in keys(QC_results)]
    SAs_all = [QC_results[name][4] for name in keys(QC_results)]
 
    return QC_results, ΣFIs_all, MFIs_all, SAs_all

end



function apply_QC(filename, results)
    # Load the QC results from the CSV file
    qc_data = CSV.read(filename, DataFrame)

    # Create an empty dictionary for storing the applied QC results
    applied_QC_results = Dict()

    for name in names(qc_data)
        # Get the QC results for the current condition
        qc_result_raw = qc_data[:, name]

        # Remove the missing values from the QC results
        qc_result = Bool.(collect(skipmissing(qc_result_raw)))

        # Get the cells for the current condition
        cells = results[name]

        # Summarize the cells with the loaded QC results
        ΣFIs, MFIs, SAs = summarize_cells(cells, qc_result)
        applied_QC_results[name] = (qc_result, ΣFIs, MFIs, SAs)
    end

    ΣFIs_all = [applied_QC_results[name][2] for name in keys(applied_QC_results)]
    MFIs_all = [applied_QC_results[name][3] for name in keys(applied_QC_results)]
    SAs_all = [applied_QC_results[name][4] for name in keys(applied_QC_results)]

    return applied_QC_results, ΣFIs_all, MFIs_all, SAs_all
end

# Plotting code
function plot_results(ΣFIs_all, MFIs_all, SAs_all)
    x = [fill(name, length(SAs_all[i])) for (i, name) in enumerate(keys(QC_results))]

    # println(x) 
    # plot ΣFIs
    pΣ = violin(x, ΣFIs_all, fillalpha=0.75, width=0.2, legend=false)
    dotplot!(x, ΣFIs_all, markeralpha=0.5, marker=(:black, stroke(0)), legend=false)
    boxplot!(x, ΣFIs_all, fillalpha=0.75, bar_width = 1, legend=false)
    xlabel!("[NIP]")
    ylabel!("Integrated tension (Atto647 ΣFI)")
    display(pΣ)
    # plot MFIs
    pM = violin(x, MFIs_all, fillalpha=0.75, width=0.2, legend=false)
    dotplot!(x, MFIs_all, markeralpha=0.5, marker=(:black, stroke(0)), legend=false)
    boxplot!(x, MFIs_all, fillalpha=0.75, bar_width = 0.2, legend=false)
    xlabel!("[NIP]")
    ylabel!("MFIs")
    display(pM)
    # plot SAs
    pS = violin(x, SAs_all, fillalpha=0.75, width=0.2, legend=false)
    dotplot!(x, SAs_all, markeralpha=0.5, marker=(:black, stroke(0)), legend=false)
    boxplot!(x, SAs_all, fillalpha=0.75, bar_width = 0.2, legend=false)
    xlabel!("[NIP]")
    ylabel!("SAs")
    display(pS)

  
end

function plot_results(QC_results)

    ΣFIs_all = [QC_results[name][2] for name in keys(QC_results)]
    MFIs_all = [QC_results[name][3] for name in keys(QC_results)]
    SAs_all = [QC_results[name][4] for name in keys(QC_results)]

    x = [fill(name, length(SAs_all[i])) for (i, name) in enumerate(keys(QC_results))]

    # plot ΣFIs
    pΣ = violin(x, ΣFIs_all, fillalpha=0.75, width=0.2, legend=false)
    dotplot!(x, ΣFIs_all, markeralpha=0.5, marker=(:black, stroke(0)), legend=false)
    boxplot!(x, ΣFIs_all, fillalpha=0.75, bar_width = 1, legend=false)
    xlabel!("[NIP]")
    ylabel!("Integrated tension (Atto647 ΣFI)")
    display(pΣ)
    # plot MFIs
    pM = violin(x, MFIs_all, fillalpha=0.75, width=0.2, legend=false)
    dotplot!(x, MFIs_all, markeralpha=0.5, marker=(:black, stroke(0)), legend=false)
    boxplot!(x, MFIs_all, fillalpha=0.75, bar_width = 0.2, legend=false)
    xlabel!("[NIP]")
    ylabel!("MFIs")
    display(pM)
    # plot SAs
    pS = violin(x, SAs_all, fillalpha=0.75, width=0.2, legend=false)
    dotplot!(x, SAs_all, markeralpha=0.5, marker=(:black, stroke(0)), legend=false)
    boxplot!(x, SAs_all, fillalpha=0.75, bar_width = 0.2, legend=false)
    xlabel!("[NIP]")
    ylabel!("SAs")
    display(pS)
  
end

function plot_results(QC_results)
    ΣFIs_all = [QC_results[name][2] for name in keys(QC_results)]
    MFIs_all = [QC_results[name][3] for name in keys(QC_results)]
    SAs_all = [QC_results[name][4] for name in keys(QC_results)]

    x = [fill(name, length(SAs_all[i])) for (i, name) in enumerate(keys(QC_results))]

    pΣ = dot_violin_boxplot(x, ΣFIs_all)
    ylabel!("Integrated tension (Atto647 ΣFI) try")
    display(pΣ)

    pM = dot_violin_boxplot(x, MFIs_all)
    ylabel!("MFIs")
    display(pM)

    pS = dot_violin_boxplot(x, SAs_all)
    ylabel!("SAs")
    display(pS)
end

function plot_results(QC_results,xs)
    ΣFIs_all = [QC_results[name][2] for name in keys(QC_results)]
    MFIs_all = [QC_results[name][3] for name in keys(QC_results)]
    SAs_all = [QC_results[name][4] for name in keys(QC_results)]

    # x = [fill(x, length(SAs_all[i])) for (i, x) in enumerate(xs)]
    x = [fill(log10(x), length(SAs_all[i])) for (i, x) in enumerate(xs)]
    @show x


    pΣ = dot_violin_boxplot(x, ΣFIs_all)
    ylabel!("Integrated tension (Atto647 ΣFI) try")
    display(pΣ)

    pM = dot_violin_boxplot(x, MFIs_all)
    ylabel!("MFIs")
    display(pM)

    pS = dot_violin_boxplot(x, SAs_all)
    ylabel!("SAs")
    display(pS)
end


function dot_violin_boxplot(x_labels, y_data)
    # p = violin(x_labels, y_data, fillalpha=0.75, width=0.2, legend=false)
    # dotplot!(x_labels, y_data, markeralpha=0.5, marker=(:black, stroke(0)), legend=false)
    # boxplot!(x_labels, y_data, fillalpha=0.75, bar_width = 0.2, legend=false, size=(1000, 600))
    p = boxplot(x_labels, y_data, fillalpha=0.75, bar_width = 0.2, legend=false, size=(1000, 600))

    xlabel!("[NIP]")
    return p
end

using CSV
using DataFrames

function save_conditions_to_csv(QC_results, filename)
    function pad_condition(condition, max_length)
        vcat(condition, fill(missing, max_length - length(condition)))
    end

    QC_df = DataFrame()
    sum_fi_df = DataFrame()
    mfi_df = DataFrame()
    sa_df = DataFrame()

    # Find the maximum length of the conditions for each metric
    max_length_QC = maximum(length.(QC_result[1] for QC_result in values(QC_results)))
    max_length_sum_fi = maximum(length.(QC_result[2] for QC_result in values(QC_results)))
    max_length_mfi = maximum(length.(QC_result[3] for QC_result in values(QC_results)))
    max_length_sa = maximum(length.(QC_result[4] for QC_result in values(QC_results)))

    for (cond_name, result) in QC_results
        QC = pad_condition(result[1], max_length_QC)
        sum_fi = pad_condition(result[2], max_length_sum_fi)
        mfi = pad_condition(result[3], max_length_mfi)
        sa = pad_condition(result[4], max_length_sa)

        QC_df[!, Symbol(cond_name)] = QC
        sum_fi_df[!, Symbol(cond_name)] = sum_fi
        mfi_df[!, Symbol(cond_name)] = mfi
        sa_df[!, Symbol(cond_name)] = sa
    end

    CSV.write(filename * "_QC.csv", QC_df)
    CSV.write(filename * "_sum_fi.csv", sum_fi_df)
    CSV.write(filename * "_mfi.csv", mfi_df)
    CSV.write(filename * "_sa.csv", sa_df)
end
