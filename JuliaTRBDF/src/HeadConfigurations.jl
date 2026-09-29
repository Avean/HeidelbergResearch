# src/HeadConfigurations.jl

# ============================================================
# Grouping of head configurations across series realizations
# ============================================================
#
# A configuration is the sorted list of head positions of one stationary
# profile, normalised to [0, 1] along the panel. Two configurations belong
# to the same group when they have the same number of heads and every head
# lies within HEAD_CONFIGURATION_TOLERANCE of the corresponding head of the
# group's reference, which is the first configuration that founded the group.
# On a periodic panel configurations that differ only by a shift are the same.

const HEAD_CONFIGURATION_TOLERANCE = 0.05
const HEAD_CONFIGURATION_LABEL_LENGTHS = 7:2:41


mutable struct HeadConfigurationGroup
    reference::Vector{Float64}
    # Normalised positions of the founding configuration; on a periodic
    # panel rotated so that its first head sits at 0.
    count::Int
end


function normalized_head_positions(
    heads::AbstractVector{DetectedHead},
    x::AbstractVector{<:Real},
    boundary_condition::Symbol,
)
    length(x) >= 2 || return Float64[]
    x0 = Float64(first(x))

    if boundary_condition == :periodic
        period = Float64(last(x) - first(x) + (x[2] - x[1]))
        return sort!([mod((head.position - x0) / period, 1.0) for head in heads])
    end

    panel_length = Float64(last(x) - first(x))
    return sort!([clamp((head.position - x0) / panel_length, 0.0, 1.0) for head in heads])
end


function _rotated_to_zero(positions::Vector{Float64}, index::Int)
    return sort!(mod.(positions .- positions[index], 1.0))
end


_circular_difference(a::Float64, b::Float64) = (d = abs(a - b); min(d, 1.0 - d))


function head_configuration_distance(
    reference::Vector{Float64},
    positions::Vector{Float64},
    periodic::Bool,
)
    length(reference) == length(positions) || return Inf
    isempty(positions) && return 0.0
    periodic || return maximum(abs.(reference .- positions))

    # The reference already has its first head at 0; try every head of the
    # candidate at 0 and keep the best match.
    return minimum(
        maximum(_circular_difference.(reference, _rotated_to_zero(positions, index)))
        for index in eachindex(positions)
    )
end


function assign_head_configuration!(
    groups::Vector{HeadConfigurationGroup},
    positions::Vector{Float64},
    periodic::Bool;
    tolerance::Float64 = HEAD_CONFIGURATION_TOLERANCE,
)
    for (index, group) in enumerate(groups)
        if head_configuration_distance(group.reference, positions, periodic) <= tolerance
            group.count += 1
            return index
        end
    end

    reference = periodic && !isempty(positions) ? _rotated_to_zero(positions, 1) : copy(positions)
    push!(groups, HeadConfigurationGroup(reference, 1))
    return length(groups)
end


function _head_configuration_cells(reference::Vector{Float64}, length_::Int, periodic::Bool)
    return periodic ?
        [mod(round(Int, p * length_), length_) + 1 for p in reference] :
        [round(Int, p * (length_ - 1)) + 1 for p in reference]
end


function _head_configuration_label(reference::Vector{Float64}, length_::Int, periodic::Bool)
    characters = fill('-', length_)

    for cell in _head_configuration_cells(reference, length_, periodic)
        characters[cell] = 'H'
    end

    return String(characters)
end


function head_configuration_labels(
    groups::Vector{HeadConfigurationGroup},
    periodic::Bool,
)
    isempty(groups) && return String[]

    # Use the shortest common length at which every group has its own label
    # and no two heads of one configuration share a character.
    for length_ in HEAD_CONFIGURATION_LABEL_LENGTHS
        labels = [_head_configuration_label(group.reference, length_, periodic) for group in groups]
        distinct = allunique(labels)
        separated = all(
            allunique(_head_configuration_cells(group.reference, length_, periodic))
            for group in groups
        )
        distinct && separated && return labels
    end

    labels = [
        _head_configuration_label(group.reference, last(HEAD_CONFIGURATION_LABEL_LENGTHS), periodic)
        for group in groups
    ]
    seen = Dict{String, Int}()

    return map(labels) do label
        seen[label] = get(seen, label, 0) + 1
        seen[label] == 1 ? label : label * " #" * string(seen[label])
    end
end


function head_configuration_order(groups::Vector{HeadConfigurationGroup})
    # Fewer heads first, then left to right.
    return sortperm(groups; by = group -> (length(group.reference), group.reference))
end
