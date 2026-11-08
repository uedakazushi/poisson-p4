# H+2Q 型：次数制限付き閉二形式の正確な目標と障害

2026-11-09 / v1。研究方針 `RESEARCH_POLICY.md` を先に読む。分類主稿の全成分数には変更なし。

## 1. 正確な複体と既存の証明範囲

H を無限遠超平面とする実退化因子 D=H+2Q（Q は既約の二次式、rank Q≥3）を考える。affine4で全多項式二形式の係数次数≤2を `Ω²_≤2` と記す。
```
Z_Q = {α∈Ω²_≤2 | Q dα-2 dQ∧α=0},
P_Q(β)=Q dβ-dQ∧β, β∈Ω¹_≤1,
E_Q=Z_Q/P_Q(Ω¹_≤1).
```
Q のdegree≤2により P_Q の像は Z_Q に入る（`d(α/Q²)=0` を満たす）。
実Poisson条件にはさらに `α²/2 = kQ² vol`, k≠0 を課す。Z_Q の全てがPoissonになるわけではない。

現行の `DOUBLE_QUADRIC_SUPPORT.md` と `RANK_THREE_QUADRIC_SUPPORT.md` の、90変数→140閉性行列に関する **保存済みQQ階数計算** は次の表を与える（今回新たに行列を再実行していない）。

| Q の射影階数と配置 | dim Z_Q | dim im P_Q | dim E_Q |
|---|---:|---:|---:|
| rank 4,5 の五配置 | 19 | 19 | 0 |
| rank 3, Hが頂点を避ける Q=u²+v²+z² | 19 | 19 | 0 |
| rank 3, 横断 Q=1+u²+v² | 20 | 19 | 1 |
| rank 3, 接触 Q=w+u² | 25 | 19 | 6 |

ゆえに「任意のrank Q≥3の閉αがdegβ≤1のP_Q(β)で表される」という無条件命題は偽。rank3の二例は**障害**であり、補集合を無視した一様消滅定理を目標にしてはならない。

## 2. 横断配置：補集合に実在するde Rham類

Q=1+u²+v² とすると α=Q du∧dv は Z_Q に入る。ω=α/Q²=du∧dv/Q は閉である。
a=u+iv, b=u-iv とおくと Q=1+ab。因子 Q=0 は C* と同型であり、ω のPoincaré留数は非零定数倍の da/a。これは非exactなので ω は A⁴\V(Q) の H²_dR に非零類を与える。従ってこの α は、次数制限を外しても rational primitive d(β/Q) による表示を持たない。

これは Z_Q/P_Q(Ω¹_≤1) の1次元の余剰について、線形代数計算とは異なる幾何学的説明を与える（余剰次元が1であること自体は既存QQ計算に依存）。
Pfaffian条件後も余剰が自動消滅するわけではない。現行 `RANK_THREE_QUADRIC_SUPPORT.md` §3 は、t=1、Ω=du∧dv+dz∧dw、β=−w dz による
α=QΩ+w dQ∧dz, α²/2=Q² vol
という実Poisson例を含む。

## 3. 接触配置：de Rhamは消えるが次数濾過が残る

Q=w+u² とすると affine complement は変数 w'=w+u² により G_m×A³ に同型。したがって H²_dR=0。
それでも E_Q は既存QQ計算で6次元。
r∈Sym²⟨v,z,w⟩（dim6）とおき、β_r=u dr-r du とすると
```
P_Q(β_r) = 2w du∧dr-u dw∧dr-r du∧dw .
```
右辺の係数次数は≤2だが、β_rの係数次数は2である。従って一般の有理de Rham exactnessでは消えない「affine-linear primitiveの次数制約」の余剰がある。これら六つがE_Qの基底となることは現行の全行列計算が与える結論であり、この式だけから独立性や全尽くしまで証明したと扱わない。

## 4. 改善の評価と次の目標

- **rank Q≥4の局所簡潔化:** 90×140行列の五配置を、一つの次数制限付き消滅定理で置換できれば実質的。さらにrank3頂点回避も同じ消滅側に入る。
- **rank Q=3の幾何学的理解:** 横断の1クラスは補集合のcohomology、接触の6クラスは濾過次数により生じる。区別しなければ新成分につながる源を誤って消す。
- **global147への影響:** H+2Qの入口は全五次因子分割の一層のみ。正規/非正規既約五次・二因子以外・rank2支持・各成分性/非包含は依然として必要。主分類証明全体がこれだけで短くなるわけではない。

数学的な次の課題は、上記 `E_Q` を極次数濾過付きde Rham/Koszul複体の上から同定し、全QのJordan標準形の選択なしに 0/1/6 を説明する**証明**を与えること。さらにPfaffian条件の切り出しと、既知成分への退化/包含の証明は別問題として扱う。
まず Q=1+u²+v² と Q=w+u² のprimitiveとresidueを直接検証し、次にrank≥4の統一消滅の必要な次数制御を調べる。計算の結論を所与の公理にして繰り返すだけの証明にしない。
de Rham類が零ならaffine-linear primitiveがあるという推論は禁止する。

## 5. 文献と既存証明

- 現行：`research/classification/DOUBLE_QUADRIC_SUPPORT.md` §§1–8、`research/classification/RANK_THREE_QUADRIC_SUPPORT.md` §§1–4。19/20/25、Pfaffian条件、退化枝の既存計算と証明を保持。
- A. Dimca, G. Sticlaru, *Koszul complexes and pole order filtrations*, Proc. Edinb. Math. Soc. 58 (2015), 333–354, arXiv:1108.3976, https://arxiv.org/abs/1108.3976 。Koszul cohomologyとpole-order filtrationを関係づける先行研究。今回の有限次数複体へ直接適用できることは未確認。
- A. Dimca, M. Saito, *Koszul complexes and spectra of projective hypersurfaces with isolated singularities*, arXiv:1212.1081, https://arxiv.org/abs/1212.1081 。今回のQ=rank3は非孤立特異点を持ちうるため、仮定の確認なしに適用しない。

今回の「留数と次数濾過という二種の障害」の説明は本ノートの導出であり、先行研究での新規性は未確認。QQ kernel次元についての今回の独立再計算は行っていない。
