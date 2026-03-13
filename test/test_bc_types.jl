@testset "Alloy BC type" begin
    ic = AlloyEquilibrium(0.8, 1.2, -0.45)
    @test ic isa AbstractInterfaceBC
    @test ic.k_partition == 0.8
    @test ic.T_m == 1.2
    @test ic.m_liquidus == -0.45

    ic2 = AlloyEquilibrium(0.7f0, 0.0f0, -0.2f0)
    @test ic2.k_partition == 0.7f0
    @test ic2.T_m == 0.0f0
    @test ic2.m_liquidus == -0.2f0
end
