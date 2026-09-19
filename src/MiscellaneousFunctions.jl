module MiscellaneousFunctions

using DimensionalData
using Statistics
import Base.Threads.@threads # For parallel computing

export nanmean, nanvar, nanstd, moving_average, bincounts


"""
    nanmean(x; dims = nothing)

Calculate the mean of array `x` while ignoring NaNs.
"""
function nanmean(x::AbstractArray; dims = nothing)
    if dims !=nothing
        return dropdims(mapslices(A -> mean(filter(!isnan, A)), x, dims=dims), dims=dims)
    else
        return mean(filter(!isnan,x))
    end
end

function nanmean(x::AbstractDimStack; dims = nothing)
    if dims !=nothing
        return maplayers(B -> dropdims(mapslices(A -> mean(filter(!isnan, A)), B, dims=dims), dims=dims), x)
    else
        return maplayers(B -> mean(filter(!isnan, B)), x)
    end
end


"""
    nanvar(x; dims = nothing)

Calculate the variance of array `x` while ignoring NaNs.
"""
function nanvar(x; dims = nothing)
    if dims !=nothing
        return dropdims(mapslices(A -> var(filter(!isnan, A)), x, dims=dims), dims=dims)
    else
        return var(filter(!isnan,x))
    end
end

function nanvar(x::AbstractDimStack; dims = nothing)
    if dims !=nothing
        return maplayers(B -> dropdims(mapslices(A -> var(filter(!isnan, A)), B, dims=dims), dims=dims), x)
    else
        return maplayers(B -> var(filter(!isnan, B)), x)
    end
end


"""
    nanstd(x; dims = nothing)

Calculate the standard deviation of array `x` while ignoring NaNs.
"""
function nanstd(x; dims = nothing)
    if dims !=nothing
        return dropdims(mapslices(A -> std(filter(!isnan, A)), x, dims=dims), dims=dims)
    else
        return std(filter(!isnan,x))
    end
end

function nanstd(x::AbstractDimStack; dims = nothing)
    if dims !=nothing
        return maplayers(B -> dropdims(mapslices(A -> std(filter(!isnan, A)), B, dims=dims), dims=dims), x)
    else
        return maplayers(B -> std(filter(!isnan, B)), x)
    end
end



"""
    moving_average(data, n; replacenans=false)

Calculate the moving average of `data` with window `n`.

Optionally replace the NaNs at the start and end with the nearest calculated value if `replacenans` is true.
"""
function moving_average(data::AbstractArray, n; replacenans=false)
    # Initialize the array
    smoothed = similar(data)
    smoothed[1:Int((n-1)/2)] .= NaN              # Put NaNs for the beginning chopped off bit
    smoothed[end-(Int((n-1)/2) - 1):end] .= NaN  # Same for the end

    smoothed[Int((n-1)/2) + 1:end-(Int((n-1)/2))] = [nanmean(@view data[i:(i+n-1)]) for i in 1:(length(data)-(n-1))]

    # optionally replace the start and end nans with the nearest calculated value 
    if replacenans
        smoothed[1:Int((n-1)/2)] .= smoothed[Int((n+1)/2)]
        smoothed[end-(Int((n-1)/2) - 1):end] .= smoothed[end-(Int((n+1)/2) - 1)]
    end

    return smoothed
end




"""
    bincounts(datax, binedgesx; makepdf=true)

Bin `datax` in 1 dimension between `binedgesx`.

Optionally normalize by bin area and number of data points to return a probability density function using `makepdf=true`.
"""
    function bincounts(datax, binedgesx; makepdf=true)

        ## The number of bins is one less than the size of the bin edges ##
        nbinsx = size(binedgesx, 1)-1

        bin_areas = diff(binedgesx)

        counts = Array{Float64}(undef, nbinsx)

        if makepdf
            @threads for x in collect(1:nbinsx)
                counts[x] = sum((datax .> binedgesx[x] .&& datax .< binedgesx[x+1])) ./ bin_areas[x]
            end
            counts = counts ./ length(datax)
        else
            @threads for x in collect(1:nbinsx)
                counts[x] = sum((datax .> binedgesx[x] .&& datax .< binedgesx[x+1]))
            end   
        end

        ## Find the bin centers for plotting ##
        bincenx = binedgesx[1:end-1] .+ diff(binedgesx)/2

        return (; counts, bincenx)
    end

"""
    bincounts(datax, datay, binedgesx, binedgesy; makepdf=true)

Impose a 2d grid on the x-y pairs of data in `datax` and `datay` and count the number of points contained in each cell.

Optionally normalize by bin area and number of data points to return a probability density function using `makepdf=true`.
"""
    function bincounts(datax, datay, binedgesx, binedgesy; makepdf=true)

        ## The number of bins is one less than the size of the bin edges ##
        nbinsx = size(binedgesx, 1)-1
        nbinsy = size(binedgesy, 1)-1

        bin_areas = diff(binedgesx) .* diff(binedgesy)'

        counts = Array{Float64}(undef, nbinsx, nbinsy)

        if makepdf
            @threads for (x, y) in collect(Iterators.product(1:nbinsx, 1:nbinsy))
                counts[x,y] = sum((datax .> binedgesx[x] .&& datax .< binedgesx[x+1]) .&& (datay .> binedgesy[y] .&& datay .< binedgesy[y+1])) ./ bin_areas[x, y]
            end
            counts = counts ./ length(datax)
        else
            @threads for (x, y) in collect(Iterators.product(1:nbinsx, 1:nbinsy))
                counts[x,y] = sum((datax .> binedgesx[x] .&& datax .< binedgesx[x+1]) .&& (datay .> binedgesy[y] .&& datay .< binedgesy[y+1]))
            end   
        end

        ## Find the bin centers for plotting ##
        bincenx = binedgesx[1:end-1] .+ diff(binedgesx)/2
        binceny = binedgesy[1:end-1] .+ diff(binedgesy)/2

        return (; counts, bincenx, binceny)
    end

end