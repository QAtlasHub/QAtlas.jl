# Identity-Based Self-Validation Harness

!!! warning "Status: Unstable (v0.18.x)"
    `verify_thermodynamic_identities` ハーネスと `SYMMETRY_IDENTITIES` 系の identity 集合は v0.17–0.18 で導入された新サーフェス。`ThermoIdentity` 構造体や `model_filter` predicate の signature は v0.19 で変更される可能性があります。テストファイル (`test/identities/`) は public API ではなく **internal verification infrastructure** であり、外部から直接 import しないでください。

## What Identities Catch

QAtlas の `fetch` メソッド群が物理的に整合しているかを **複数の独立な経路で同じ量を計算して一致するか** で機械的に検証。観測量実装の sign error / convention drift / kwarg drift を ED ベースの単純比較より早く / 広く検出。

## 6-axis Test Diversification

QAtlas v0.18 では以下 6 軸の identity 検証層が並走する:

| 軸 | 例 | 関連 PR |
|---|---|---|
| (A) Cross-method | `c_v` を `Var(H)` / `∂ε/∂β` / `∂s/∂β` の 3 経路で算出して一致 | #132 |
| (B) Cross-symmetry | SU(2) 不変 model で `χ_xx = χ_yy = χ_zz`, `m_α = 0` | #133 |
| (C) Cross-model | TFIM `J=0` / `h=0` の textbook closed form 帰着 | #133 |
| (D) Cross-BC scaling | OBC/PBC `→` Infinite の `1/N`, `1/N²` 指数を fit | #135 |
| (E) Property-based | random sweep + 不変量 (`Var ≥ 0`, `s ≥ 0`, `0 ≤ S_α ≤ ℓ log 2`) | #136 |
| (F) Time-domain | dynamic 相関子の reflection symmetry / Hermiticity / locality | #134 |

## Core API: `ThermoIdentity` and `verify_thermodynamic_identities`

(`test/util/thermodynamic_identities.jl` 参照)

```julia
struct ThermoIdentity
    name::String
    requires::Vector{Type}              # required quantity types for dispatch lookup
    check::Function                     # (model, bc, params) -> (lhs, rhs)
    model_filter::Function              # default: m -> true
end

verify_thermodynamic_identities(model, bc; βs, identities=DEFAULT_IDENTITIES, rtol=1e-8, atol=1e-10)
    -> Vector{IdentityCheckResult}
```

`IdentityCheckResult.status` は `:pass`, `:fail`, `:skipped` のいずれか。`:skipped` は `requires` の dispatch が無いか `model_filter` が false を返したケース (NaN, NaN)。

## DEFAULT_IDENTITIES (4 universal)

1. `GIBBS_RELATION` — `ε = f + T·s`
2. `SPECIFIC_HEAT_FROM_ENERGY` — `c_v = -β² ∂ε/∂β` (ForwardDiff)
3. `SPECIFIC_HEAT_FROM_ENTROPY` — `c_v = -β · ∂s/∂β` (ForwardDiff)
4. `MAGNETIZATION_X_FROM_FREE_ENERGY` — `m_x = -∂f/∂h` (central diff; skipped if model has no `h` field)

## SYMMETRY_IDENTITIES (5 opt-in)

- `SU2_CHI_XX_EQ_YY`, `SU2_CHI_YY_EQ_ZZ` — SU(2) 不変点で χ 軸等値
- `MAGNETIZATION_X_VANISHES_SU2`, `MAGNETIZATION_Z_VANISHES_SU2` — m_α = 0 (canonical, SU(2))
- `MAGNETIZATION_Y_VANISHES_REAL_H` — 実 Hamiltonian で m_y = 0 (parity)

`is_su2_symmetric(model)` predicate でモデル ごとに有効化:
- `Heisenberg1D` → true
- `S1Heisenberg1D` → true
- `XXZ1D` → `isapprox(Δ, 1.0)`
- 他 → false

## SUSCEPTIBILITY_XX_KUBO_FROM_MAGNETIZATION (opt-in only)

`χ_xx = ∂m_x/∂h` の Kubo 静的応答テスト。**OBC dense-ED 系 backend は equal-time variance 規約で `β·Var(M_x)/N` を返す**ため Kubo identity と不一致 — DEFAULT には入れない。Calabrese-Mussardo 系の closed-form (TFIM Infinite) でのみ pass。詳細は `test/util/thermodynamic_identities.jl` の docstring と PR #132 のレビュー。

## Kubo vs Variance — Convention Note

OBC `χ_xx` (TFIM/XXZ1D/S1Heisenberg) は equal-time variance:

    χ_xx^(var) = β · Var(M_x) / N = (β/N) · Σᵢⱼ [⟨σˣᵢσˣⱼ⟩ − ⟨σˣᵢ⟩⟨σˣⱼ⟩]

Infinite TFIM `χ_xx` は Calabrese-Mussardo Kubo 形:

    χ_xx^(Kubo) = ∂m_x/∂h = (1/Nβ) ∂²(log Z)/∂h²

古典極限 `[M_x, H] = 0` では一致するが、量子では operator-ordering 補正で異なる。harness はこの違いを検出可能 (PR #132)。

## Per-model Identity Test Files

各モデルの test/identities/test_identities_*.jl で harness を呼ぶ:

- `test_identities_TFIM.jl`, `test_identities_TFIM_pbc.jl`
- `test_identities_XXZ1D.jl`
- `test_identities_S1Heisenberg1D.jl`
- `test_identities_Heisenberg1D.jl`
- `test_identities_KitaevHoneycomb.jl`
- `test_TFIM_dynamic_symmetries.jl` — (F) time-domain
- `test_TFIM_limits_cross_model.jl` — (C) cross-model
- `test_cross_bc_scaling.jl` — (D) BC scaling
- `test_property_invariants.jl` — (E) property-based

合計 **約 1000+ identity-style tests** が CI で走る。

## Adding a New Identity

```julia
const MY_IDENTITY = ThermoIdentity(
    "lhs = rhs description",
    Type[Quantity1, Quantity2],
    function (model, bc, params)
        β = params.β
        lhs = fetch(model, ...)
        rhs = fetch(model, ...)
        return Float64(lhs), Float64(rhs)
    end;
    model_filter=is_su2_symmetric,  # optional, default _ -> true
)
```

then either add to `DEFAULT_IDENTITIES` (universal) or pass via `identities=[MY_IDENTITY]` kwarg.

## References

- PR #126 (initial harness, #117 issue)
- PR #132 (cross-method, Kubo vs variance discovery)
- PR #133 (symmetry + cross-model)
- PR #134 (time-domain identities)
- PR #135 (cross-BC scaling)
- PR #136 (property-based invariants)
