# models/quantum/TFIM/TFIM_random_registry.jl
#
# Random transverse-field Ising chain — the exact Griffiths dynamical exponent
# and the fixed point it flows to at J == h.

@register(
    RandomTFIM,
    DynamicalExponent,
    Infinite,
    method=:analytic_griffiths_root,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/TFIM/test_TFIM_random_griffiths.jl",
    references=["IgloiMonthus2005", "FisherDS1995"],
    notes="Positive root of [(J/h)^{1/z}]_av = 1 (Eq. 4.15); on this disorder family (J/h)^{1/z} = 1 - D^2/z^2, z > D. Throws at J == h, where no finite z exists.",
)

@register(
    RandomTFIM,
    ActivatedExponent,
    Infinite,
    method=:analytic,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/TFIM/test_TFIM_random_griffiths.jl",
    references=["IgloiMonthus2005", "FisherDS1995"],
    notes="psi = 1/2 at the infinite-randomness fixed point (Eq. 4.13); refused off criticality, where a finite z is the operative exponent instead.",
)

@register(
    RandomTFIM,
    UniversalityClass,
    Infinite,
    method=:analytic,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/TFIM/test_TFIM_random_griffiths.jl",
    references=["IgloiMonthus2005", "FisherDS1995", "RefaelMoore2004"],
    notes="Universality class identity at J == h: :IsingSDRG (IRFP).",
)
