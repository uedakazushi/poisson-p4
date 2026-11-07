# 非正規二重曲面のglobal degree制限とrank二CI(2,2)枝の排除

2026-11-07、主義務G。DをP4の既約・被約五次とし、Sをその非正規曲面の既約成分とする。
第1節はPoisson rankを仮定しない。
第2–4節ではSが滑らかで一般ordinary node、π|Sが一般rank二とするが、
normalization X又はconductorのsmoothnessは仮定しない。

## 1. hyperplane内のdegree三以上の曲面は全て不可能

Sがhyperplane内のdegree e≥3のintegral surfaceなら
`I_S=(ℓ,q_e)` はprime complete-intersection idealである。
DのS一般点での重複度は二以上なので、定義五次式Fは `I_S^(2)` に入る。
prime complete intersectionのsquareはprimaryであり、従って `I_S^(2)=I_S²`。
これは例えば生成元(ℓ²,ℓq,q²)のHilbert--Burch行列

`[−q,0; ℓ,−q; 0,ℓ]`

によるgrade二のperfect性から分かる。quotientがCMでembedded primeを持たず、
唯一のminimal primeはI_SなのでI_S-primaryとなる。

2e>5のため、degree五では

`F=ℓ²A₃+ℓq_e B_{4−e}`。

従ってℓがFを割り、既約性に矛盾する。
Sのnormality、smoothness、πのrankは不要で、全退化したq_eもintegralである限り含む。
よってhyperplane内の非正規曲面としてはplane又はquadricだけが残る。

## 2. rank二nodeのhorizontal/vertical分解と実anticanonical曲線

Sが滑らか、DがS一般点でordinary node、π|Sが一般rank二とする。
`E=P(N_S/P4)`、`ξ=c₁(O_E(1))`、H=O_S(1)とする。
`IRREDUCIBLE_NONNORMAL_RANK2_CLEAN_SURFACE.md` のblowup/residue議論は、
vertical成分がある場合もhorizontal部分にそのまま適用できる。
E上の閉residue一形式のhorizontal biresidueは非零定数λと−λなので、
horizontal部分は相異なる二degree一divisor Γ₁,Γ₂にsplitする。

全normal-cone divisorを

`D'|E=Γ₁+Γ₂+p*V`

と書く。VはS上のeffective curve divisorで、multiplicityも保持する。
`O_E(Γᵢ)=O_E(1)⊗p*Aᵢ` とすると

`A₁+A₂+V=5H`。

Γᵢはnormal線形form `ℓᵢ:N→Aᵢ` で定義される。
vertical factorを最大限分離したので、ℓᵢの係数には共通divisorial factorが無い。
normal-cone五次sectionは `γℓ₁ℓ₂`、div γ=V。
ordinary nodeのため `δ=det(ℓ₁,ℓ₂)` は零でない。
δは `A₁+A₂−det N=−K_S−V` のsectionである。

π|Sのanticanonical sectionをσとする。
normal座標を(ℓ₁,ℓ₂)へ変えるとPfaffianのnormal Jacobianはδ倍となる。
log skew行列の逆で接方向Poisson成分は、biresidue λとlog Pfaffianの逆の積である。
従ってorientationによる符号を除き

`σ=λ γ δ`。

これはgeneric ordinary-node開上の恒等式であり、両辺が同じline bundleの大域sectionなので全Sへ延びる。
従って実zero divisor R_πについて

`R_π=V+div δ`、特に `V≤R_π`。                          (1)

これはclassだけの一致ではない。
verticalでのpoleが一位とは仮定せず、horizontal一般点の非零定数biresidueだけを使った。

## 3. smooth quartic del Pezzoのvertical全枝は可約quinticになる

S=V(q₁,q₂)⊂P4がsmooth complete intersection of two quadricsなら
`−K_S=H`, `H²=4`, `N=O_S(2H)²`。
第2節のnormal線形formの係数を

`ℓᵢ=(sᵢ,tᵢ)`、`Bᵢ=Aᵢ−2H`

とするとsᵢ,tᵢはBᵢの二sectionで、共通fixed divisorを持たない。
従ってBᵢ=0又はBᵢは非零mobile nef classである。
また

`B₁+B₂+V=H`。                                           (2)

非零mobile BについてHB≥2である。
実際HB=1ならHodge indexから `0≤B²≤1/4` なのでB²=0。
Riemann--Rochのparity `B²−K_SB=B²+HB` が偶数であることに矛盾する。

V≠0の場合、両Bᵢが非零なら(2)へHを掛けて
`4=HB₁+HB₂+HV≥2+2+1` となって矛盾する。
従って一方のBᵢは零。そのℓᵢは二normal quadricのconstant linear combinationである。
constant GL₂変更によりそのfactorをq₂とする。

prime CIのsquare/symbolic-square一致とdegree五から

`F=q₁²L₁+q₁q₂L₂+q₂²L₃`

であり、各Lᵢはambient linear formである。
normal-cone先頭がq₂で割れるので `L₁|S=0`。
SはP4で非退化なのでL₁=0、従ってF全体がq₂で割れる。
I_S³はdegree六未満の元を持たず、ここに未記録の高normal-order項を加える余地は無い。
既約quinticに矛盾する。

## 4. 同じCI(2,2)のvertical無し全枝も不可能

V=0では閉residue一形式のGysin関係から `A₁≡A₂`。
和は5Hなので、整数line bundle A=A₁−2Hについて `2A≡H`、A²=1。
K_S=−HなのでK_SA=−2。
するとRiemann--Rochの `A²−K_SA=3` は奇数となって矛盾する。

従って**smooth CI(2,2)を一般rank二ordinary-double曲面に持つ既約五次Poisson退化因子は存在しない**。
vertical有無の全補集合を覆い、normalization/conductorのglobal smoothnessは使わない。
rank零、非node、特異CI又は他の曲面についてこの命題を無断で拡張しない。

## 5. 次の有限入口

S=planeのvertical枝V≠0ではN=O(1)²、E=P2×P1で、horizontal Γᵢのbidegreeは `(aᵢ−1,1)`、aᵢ≥1。
`v+a₁+a₂=5` と(1)から `1≤v≤3`。
v=1,2ではσのcubicはvertical factorと非定数det factorに分かれる。
v=3ではhorizontal二方向はconstantで、σの全cubicがverticalになる。
これは全vertical plane枝の有限次数入口であり、そのPoisson源全帰属又は非優越まで証明したものではない。
quadric、他のrational曲面、normal/non-normal特異曲面も別に残す。
