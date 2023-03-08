
"""
A dictionary that derives from another one but with (periodically)
shifted indices.
"""
struct IndexShiftedDict{S,T} <: Dictionary{S,T}
    superdict   ::  Dictionary{S,T}
    shift       ::  Int
end

superdict(Φ::IndexShiftedDict) = Φ.superdict

Base.size(Φ::IndexShiftedDict) = size(superdict(Φ))

BasisFunctions.support(Φ::IndexShiftedDict) = support(superdict(Φ))

shift_index(Φ::Dictionary, idx, shift) =
    shift_index(eachindex(Φ), idx, shift)
shift_index_back(Φ::Dictionary, idx, shift) =
    shift_index_back(eachindex(Φ), idx, shift)

shift_index(Φ::IndexShiftedDict, idx) =
    shift_index(eachindex(Φ), idx, Φ.shift)
shift_index_back(Φ::IndexShiftedDict, idx) =
    shift_index_back(eachindex(Φ), idx, Φ.shift)

function shift_index(range::AbstractUnitRange, idx, shift)
    idx2 = idx - shift
    if idx2 < first(range)
        idx2 + length(range)
    else
        idx2
    end
end

function shift_index_back(range::AbstractUnitRange, idx, shift)
    idx2 = idx + shift
    if idx2 > last(range)
        idx2 - length(range)
    else
        idx2
    end
end

import BasisFunctions:
    basisfunction,
    unsafe_eval_element,
    unsafe_eval_element_derivative

basisfunction(dict::IndexShiftedDict, idx) =
    basisfunction(superdict(dict), shift_index_back(dict, idx))

unsafe_eval_element(dict::IndexShiftedDict, idx, x) =
    unsafe_eval_element(superdict(dict), shift_index_back(dict, idx), x)
unsafe_eval_element_derivative(dict::IndexShiftedDict, idx, x, order) =
    unsafe_eval_element_derivative(superdict(dict), shift_index_back(dict, idx), x, order)
