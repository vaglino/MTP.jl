"""
Calculates the mean fluorescence intensity (MFI) of tension of a Cell
# Arguments
- `cell`: Cell object
- `channel`: channel to be used for MFI calculation (default = 1 for MTP, 2 for RICM)
# Returns
- `MFI`: mean fluorescence intensity of tension
"""
function MFI(cell::Cell; channel = 1)
    F = cell.I[:, :, channel]
    seg_mask = cell.seg[:, :, 2] .!= 0
    # imshow(F .* seg_mask)
    if spread_area(cell) != 0.0
        MFI = sum(F .* seg_mask) / spread_area(cell)
    else
        MFI = 0.0
        # seg_mask = cell.seg[:,:,3] .!= 0
        # MFI = sum(F .* seg_mask) / spread_area(cell)
    end
    # MFI = mean(F[seg_mask])

end

"""
Calculates the sum of fluorescence intensity (ΣFI) of tension of a Cell
# Arguments
- `cell`: Cell object
- `channel`: channel to be used for ΣFI calculation (default = 1 for MTP, 2 for RICM)
# Returns
- `ΣFI`: sum of fluorescence intensity of tension
"""
function ΣFI(cell::Cell; channel = 1)
    F = cell.I[:, :, channel]
    seg_mask = cell.seg[:, :, 2] .!= 0
    # imshow(F .* seg_mask)
    ΣFI = sum(F .* seg_mask)
end

"""
Calculates the spread area of a Cell
# Arguments
- `cell`: Cell object
- `channel`: channel to be used for spread area calculation (default = 1 for MTP, 2 for RICM)
# Returns
- `SA`: spread area of cell
"""
function spread_area(cell::Cell)
    # F = cell.I[:,:,channel]
    seg_mask = cell.seg[:, :, 2] .!= 0
    seg_mask = morph.binary_erosion(seg_mask, morph.disk(5))
    # imshow(seg_mask)
    SA = sum(seg_mask) * px_area
end

"""
Calculate MFI, ΣFI and SA for all cells in a given set of cell images
# Arguments
- `cells`: array of Cell objects
- `QC`: array of booleans indicating which cells to include in analysis
# Returns
- `ΣFIs`: array of ΣFI values
- `MFIs`: array of MFI values
- `SAs`: array of SA values
"""
function summarize_cells(cells, QC)
    cells_QC = cells[Bool.(QC)]
    ΣFIs = ΣFI.(cells_QC; channel = 1)
    MFIs = MFI.(cells_QC; channel = 1)
    SAs = spread_area.(cells_QC)
    return ΣFIs, MFIs, SAs
end

"""
Calculates the area of a Cell based on DIC segmentation
# Arguments
- `cell`: Cell object
- `channel`: channel to be used for area calculation (default = 1 for MTP, 2 for RICM)
# Returns
- `cell_area`: area of cell
"""
function cell_area(cell::Cell)
    # F = cell.I[:,:,channel]
    seg_mask = cell.seg[:, :, 1] .!= 0
    # imshow(seg_mask)
    cell_area = sum(seg_mask) * px_area
end
