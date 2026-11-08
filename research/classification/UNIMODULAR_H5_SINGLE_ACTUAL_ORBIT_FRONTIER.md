# 全unimodular h5支持：既知閉和と一つの実軌道への有限被覆

2026-11-08、PRIMARY `E_multiple_hyperplane`。producer純合成、独立/root査読前。
新ambient/operator/minor/point/jet run無し。

## 1. 主張と正確な残義務

全quadratic affine Poisson tensor P=C+L+Q、divP=0、Pf(P)=nonzero constantについて、
斉次化した全projective点は、有限個の既知閉成分の和と、次の**一つのactual PGL軌道**の和に入る：

```
P_*=(x,−x²,xz−1,xy−1,−yz+w−y,−xw+xy+y+z),
pair order=(xy,xz,xw,yz,yw,zw), Pf=1, div=0.
```

この主張はP_*の既知親又は成分性をまだ決めない。
従って全h5 Eの完成を宣言せず、その一軌道のmembership／全ambient親を唯一の残義務として残す。
Gの既完了非優越を再開しない。全scheme厚み／B18は対象外。

## 2. 全最高部と全低次fiberの入口

volume dualな閉二形式α=J+β+γで、全C6/L24/Q60と全constant-Pf条件を使う。
θ=i_Eγ/4、dθ=γの線形同定により、全最高部はPym Theorem3.1と直後のLPT拡張の
六葉層閉包L1111、L112、R22、R13、S23、E3に入る。
この外部入力は既にGで独立scope確認済みであり、今回の純証明で再証明したとはしない。
θの共通factorを持つ補集合、γ零、γのfactor付きも保持する。

各全閉包に対する適用は次の通り。

| 全最高部支持 | 全低次延長の現在の矢印 |
|---|---|
| γ零 | primitive isotropic βからconstant Hamiltonian、全CA |
| γ factor付き | 線形ηのGrassmannian二rulingからlinear-Casimir又はfullpullback、双方既知被覆 |
| R13全閉包 | proper linear-Casimir incidence、全affine Hamiltonian対とliteral modular化曲線 |
| R22全閉包 | 実Grassmannian pencil全像。primitive regularは全T16/空、factorは上行、全退化pencilのcommon-kernelはfullpullback、残Kroneckerは全空 |
| S23全閉包 | proper constant-vector horizontal incidence、fullpullbackの全C/L定理 |
| L1111/L112全閉包 | proper quartic integrating-factor incidenceと全pole型、全sourceは既知又はP_*一軌道 |
| E3全閉包 | proper nonzero nilpotent tangent-incidence全四Jordan型。主軌道空、regular境界はlinear-Casimir、rank二はclosedL112、rank一はpullback、零も上行 |

R22、S23、E3の閉包の橋と全nil kernelは
`REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` §16、
`UNIMODULAR_EXCEPTIONAL_NILPOTENT_INCIDENCE.md` の既全source定理を使う。
新しいE結論は、各点における**全C/L extension**に対して以下の全membershipを適用した結果である。

* `UNIMODULAR_LINEAR_HAMILTONIAN_WHOLE_MEMBERSHIP.md`。
* `UNIMODULAR_REGULAR_PENCIL_WHOLE_T16_MEMBERSHIP.md`。
* `UNIMODULAR_ZERO_AND_FACTOR_TOP_ENTRANCE.md`。
* `UNIMODULAR_PULLBACK_WHOLE_KNOWN_COVERAGE.md`。
* `UNIMODULAR_LOG_POLE_WHOLE_BETA_ZERO_AND_SINGLE_ORBIT.md`。

最高Qのgeneric平滑化から未知fiberの被覆を推測せず、
各下位証明は全低次源からのliteral curve、全必要十分linear source、又は全unit-Pf残式を使う。

## 3. 残軌道のsource完全性

P_*はx²y² integrating-factor境界のaligned-zero源で、旧sourceの全残式とPf非零は
n0=n3=a0=0、n1*n2≠0を強制する。
全四translationと最高GL4/tensor normalizerを戻し、
explicit diagonal map/tensor scaleでn1=n2=1へ移せることを新log-pole noteで証明した。
これは一つのpointを試して全sourceとみなしたのではなく、
全sourceのactual正逆を先に証明してforced canonical pointを得たものである。

P_*の全ambient判定を必要とする場合も、この同じ一点だけを使う。
別point、旧Gのrank不足の再試験、adaptive jet/minor列へ進む理由は無い。
局所判定を使うなら事前の有限contractと、smooth-positive／非閉鎖のscope分離を要する。
既知閉親へのactual曲線又は全係数同定でこの軌道を覆えれば、新ambient算術自体は不要である。
