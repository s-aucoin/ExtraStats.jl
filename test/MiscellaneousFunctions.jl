@testset "MiscellaneousFunctions" begin
    mod = ExtraStats.MiscellaneousFunctions

    v = [1.0, NaN, 3.0, 5.0]
    @test mod.nanmean(v) ≈ 3.0
    @test mod.nanvar(v) ≈ 4.0
    @test mod.nanstd(v) ≈ 2.0

    smooth = mod.moving_average([1.0, 2.0, 3.0, 4.0, 5.0], 3)
    @test isnan(smooth[1])
    @test smooth[2:4] ≈ [2.0, 3.0, 4.0]
    @test isnan(smooth[end])

    counts = mod.bincounts([0.1, 0.2, 1.7, 1.9], [0.0, 1.0, 2.0]; makepdf=false)
    @test sum(counts.counts) == 4

    counts2 = mod.bincounts([0.1, 0.2, 1.7, 1.9], [0.1, 0.9, 1.8, 2.0], [0.0, 1.0, 2.0], [0.0, 1.0, 2.0]; makepdf=false)
    @test size(counts2.counts) == (2, 2)
end
