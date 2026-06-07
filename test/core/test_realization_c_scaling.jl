# Systematic (rough) verification of the model ↔ universality correspondence:
# at each registered critical `example`, the model's OWN entanglement-entropy
# scaling must track the realized class's central charge. We compare the
# *difference* S(ℓ₂) − S(ℓ₁) of the model's EE against the universal Calabrese–
# Cardy formula evaluated with the class c (the non-universal constant cancels).
# Loose tolerance — this is a "does it scale with the right c" sanity sweep, not
# a precision c-extraction (that is the deferred verification-density work).

using QAtlas, Test
using QAtlas: REALIZES, Universality, realized_class

@testset "realization c-scaling (systematic, rough)" begin
    ℓ1, ℓ2 = 20, 80   # subsystem sizes (Int — model EE fetches expect integer ℓ)
    checked = 0
    for r in REALIZES
        r.example === nothing && continue
        u = Universality(r.class)
        # class must be a 1+1D CFT (central charge defined), else skip
        c = try
            QAtlas.fetch(u, CentralCharge())
        catch
            continue
        end
        # the model must compute its own infinite-system VonNeumannEntropy, else skip
        local sm1, sm2, su1, su2
        try
            sm1 = QAtlas.fetch(r.example, VonNeumannEntropy(), Infinite(); ℓ=ℓ1)
            sm2 = QAtlas.fetch(r.example, VonNeumannEntropy(), Infinite(); ℓ=ℓ2)
        catch
            continue
        end
        su1 = QAtlas.fetch(u, VonNeumannEntropy(), Infinite(); ℓ=ℓ1)
        su2 = QAtlas.fetch(u, VonNeumannEntropy(), Infinite(); ℓ=ℓ2)
        Δm, Δu = sm2 - sm1, su2 - su1
        checked += 1
        # the model's critical EE should scale with the class c (same log slope)
        @test isapprox(Δm, Δu; rtol=0.25)
        # and its own predicate resolves back to this class
        @test realized_class(r.example) === r.class
    end
    @info "c-scaling verified on $(checked) realization example(s)"
    @test checked >= 1   # at least one CFT realization is systematically checked
end
