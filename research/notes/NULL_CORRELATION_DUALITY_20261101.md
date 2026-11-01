# Null-correlation 曲線、双対楕円五次曲線、および無限遠の同変性

2026-11-01。R0 型無限遠の独立監査として、primitive 第一ジェット β に付随する
P³ 内の零曲線をベクトル束として同定する。さらに、既存の Cβ⊂LG(2,4) と
**直交する Grassmannian 切断**であることを全16変数の恒等式で示す。
この曲線は滑らかな範囲で Pic³(Cβ) であり、Jacobian は Cβ と同型である。
一般関数体上で Cβ 自体と同型であるとは主張しない。

## 1. 16次元の section 空間

V を4次元、J を標準 Poisson 行列、Ω=J⁻¹=−J を逆2-form とする。
P(V) 上の null-correlation bundle N は

    0 → O(−1) → Ω¹_P³(1) → N → 0

で定義する。左の写像は Ω によるもの。2だけ捻ると

    0 → O(1) → Ω¹_P³(3) → N(2) → 0.                 (1)

primitive な閉線形2-form β に対し、α=ι_Eβ=−β(y)y と置く。
α は係数が二次で Euler vector E による contraction が零なので、Ω¹(3) の section である。
その N(2) への像 τβ の零点条件は

    τβ([v]) = β_v(v,−)|_{v⊥/⟨v⟩} = 0.               (2)

これは α(v) が Ω(v,−) に比例する条件であり、quotient (1) の零条件と一致する。

この線形写像 W16→H⁰(N(2)) は同型である。
Euler sequence より h⁰(Ω¹(3))=4·10−20=20、(1) より h⁰(N(2))=20−4=16。
一方、α の類が零なら α=fι_EΩ、f は線形である。Cartan の公式から

    3β=dα=df∧ι_EΩ+2fΩ.

symplectic trace を取ると 5f/3 となるため、primitive 条件は f=0 を強いる。
従って β=0 であり、写像は単射。同次元なので同型となる。
コードではこの Cartan 恒等式と、任意の線形 f に対する trace を記号的に検証した。
固定行列で成分 β01+β23 を取る規約では trace は −5f/3 となる。

## 2. 一般零曲線の次数・種数・連結性

Ω¹(3) は大域生成されるので N(2) も大域生成される。
c₁(N)=0,c₂(N)=H² より

    c₁(N(2))=4H,    c₂(N(2))=5H².

標数零で一般 section の零集合は滑らかな余次元2の曲線で、次数5を持つ。
空集合でないことは c₂≠0 からも従う。adjunction は ω_E≅O_E を与える。
連結性も必要なので、Koszul sequence

    0 → O(−4) → N(−2) → I_E → 0

を用いる。Euler sequence と null-correlation sequence により N(−2) は acyclic。
また H²(O(−4))=0 なので H¹(I_E)=0、従って h⁰(O_E)=1。
以上で一般零曲線 Eβ は滑らかな連結な種数1・次数5の曲線となる。
これは P³ 内では linearly normal でない。実際 h⁰(O_E(1))=5 だが、周囲の
H⁰(O_P³(1)) は4次元である。

## 3. 明示的な第五座標と Pfaffian model

第五座標を t として P⁴ 内の交代行列を

    Φ_E(y,t) = [ β(y)+tJ    Jy ]
               [ −(Jy)ᵀ     0 ]                         (3)

と置く。符号付き4次 Pfaffian ベクトルは正確に

    ( −J(β(y)−tJ)y,  t²+Pf(β(y)) ).                     (4)

primitive な4次交代行列の恒等式 (Jβ)²=−Pf(β)I から、最初の四式を満たす点では

    (t²+Pfβ)y=0

が従う。従って y のいずれかの座標が可逆な各 chart では第五式も従う。
y=0 では第五式が t²=0 になるため、(3) の Pfaffian locus は射影中心 [0:0:0:0:1] を含まない。

射影 [y:t]↦[y] は (2) の零スキームへの同型である。
実際、(Jy)_i≠0 の chart で

    t=(β(y)y)_i/(Jy)_i

が逆写像を与え、四式 (β−tJ)y=0 は α と Jy の比例条件に一致する。
第五式は上の行列恒等式で従う。この chart ごとの議論は被約点だけでなくスキームとして成立し、
β の一般性を必要としない。

## 4. Cβ との直交双対性

既存の Cβ の Grassmannian 表示で、primitive な bivector Z∈Λ²V、
ℓβ(Z)_k=Σ_{i<j}β_{ij,k}Z_{ij} に対して

    Kβ(Z) = [ Z        −Jℓβ(Z) ]
            [ (Jℓβ(Z))ᵀ    0   ]

とする。Z は ⟨J,Z⟩=0 を満たし、この graph は Λ²(V⊕k) 内の5次元部分空間である。
(3) の係数行列も5次元部分空間を張る。全順序なし添字の自然な pairing について

    ⟨Kβ(Z),Φ_E(y,t)⟩
      = ℓβ(Z)·y + (Jy)ᵀ(−Jℓβ(Z)) = 0.                 (5)

従って (3) の五次元空間は正確に Kβ の annihilator である。
既存の canonical matrix Φ_C=[[star Z,ℓ],[-ℓᵀ,0]] との規約も

    diag(J,1) Φ_C diag(J,1)ᵀ = −Kβ(Z)

で一致する。(4)、(5)、この congruence は β の全16変数を残して検証した。

ここで一次文献 Shinder–Zhang の Proposition 2.8 を適用する。同定した直交線形切断の
一方 Cβ が滑らかなら他方も滑らかで、後者は Pic³(Cβ) となる。
したがって Δ≠0 の範囲で Eβ≅Pic³(Cβ)。次数5の偏極を用いた双対化により Pic² とも同型。
Jacobian は Cβ の Jacobian と同型である。複素数上の各滑らかな fiber は抽象曲線として
Cβ と同型だが、一般関数体では異なる torsor になり得るため同一視しない。
この文献上の双対性は本計算の新規性の主張ではない。
[Shinder–Zhang, Proposition 2.8](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/blms.12339).

二つの有理第一ジェットで Fisher の c4,c6 が既存 Φ_C と正確に一致することも再計算した。
これは符号・正規化の点検であり、普遍的双対性の証明は (5) と引用定理による。

## 5. G6 を含む135式の leading locus の同変性

`S61_UNIVERSAL_WEIGHT6_REVIEW_20261031.md` §6 の証明を独立に読み直した。
そこで flag に適用する前に証明されている事実は、より一般に次の通りである。

    g·G6_top − G6_top ∈ (旧134式の q-最高二次部分),  g∈Sp₄.

理由は、70本の最高式の q-linear syzygy 空間 S の次元が36、通常の Bianchi 部分 B が
35次元であり、S/B が Sp₄ の自明な一次元表現となることである。
同変な低次 Jacobi の特解 Q₀ を選び、全 Bianchi 恒等式

    [q,J₄]=−[Q₀(β),J₄]−[L(β),J₅]

を使うと変換差は q に依存しない係数による β²F₄+βF₅ 型の結合となる。
primitive gauge と旧 gauge の差も同じ型で処理される。
この議論は q が flag 型であることを使わない。
従って **135本の最高二次式の零スキーム全体**が Sp₄ に安定であり、R0 の非flag族にも
同じ移送が使える。これは既存の普遍 syzygy 証明の適用範囲の点検である。
今回のプログラムが36次元核の計算を新たに独立実行した、という意味ではない。

この結果は未飽和斉次化の無限遠に関するものである。t飽和された真の射影閉包、
アフィン枝の存在、全一般ファイバーの有限性・尽くしは含まない。

```sh
.venv/bin/python code/null_correlation_duality_20261101.py \
  --output /tmp/null-correlation-duality-review.json
```

証明書：`reports/NULL_CORRELATION_DUALITY_20261101.json`。
束の cohomology・Bertini・adjunction の証明は本ノート、普遍的行列恒等式はコード、
Picard 同定は原典 Proposition 2.8 と、それに必要な直交性の恒等式が根拠である。
