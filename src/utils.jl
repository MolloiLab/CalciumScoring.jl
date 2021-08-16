"""
	mask_elements(img, threshold, cutoff)
Remove every element that is below the `threshold` and
adjacent to more than the `connectivity` amount of other 
elements similarly below the threshold. Returns a mask where 
all elements are 1 for elements that should be kept based on
`threshold` and `connectivity` criteria.

# Arguments
- img: array which will be used to determine mask
- threshold: cutoff value for elements that should be 
    removed
- connectivity: minimum number of adjacent elements that are 
    required to remove an element
"""
function mask_elements(img, threshold, connectivity)
    mask = img .< threshold
    label = Images.label_components(mask)
    n = length(unique(label))
    indices = []
    for lbl in 0:(n - 1)
        lbl_mask = map(x -> x == lbl, label)
        num = count(lbl_mask)
        if num < connectivity
            rmv_indices = findall(x -> x == lbl, label)
            push!(indices, rmv_indices)
        end
    end
    for i in 1:length(indices)
        idx = indices[i]
        mask[idx] .= 0
    end
    return iszero.(mask)
end
