# 大域構造研究ノート

2026-11-09 / v3。最初に RESEARCH_POLICY.md、次に GLOBAL_PLAN.md を読む。

## 1. 現行の証明とその場所

原稿 `papers/poisson_divisor_symmetries_ja.tex` は独立に読める。詳細証明をここに重複しない。

| 主張 | 仮定・結論 | 状態／所在 |
|---|---|---|
| 因子由来可換対称性 | P^(2n)、generic nondegenerate、Dredがr因子 ⇒ r−1次元の可換対称性。非被約可 | 証明収録・今回再読、thm:abelian。新規性未確認 |
| 軌道上界 | P4でorbit≤25−r。等号・成分性は別 | cor:orbit |
| 固定因子復元 | d_Fα=0 と qα∈C×F、actual multiplicityを保持 | thm:fibre |
| 全rank4入口 | Hあり／HなしQ+C／既約五次 | thm:three-entrances |
| 標識超平面 | h-central lift一意、affine二次模型。mult_H D=5−deg f | thm:marked-affine |
| modular復元 | 固定f≠0,Zでdiv P=Z、P#df=fZ、Pf P=f iff全Poisson | thm:modular-reconstruct。記号恒等式も確認 |
| 共鳴 | Z=S+Nから[S,P]=[N,P]=0、有限な単項式共鳴型 | prop:resonance。有限な固有値比とは言わない |
| 相対復元・持上げ | 固定Fのkernelに加え、因子も動く一次・二次式 | prop:relative, prop:lifting。全高次無障害は未主張 |
| グラフ閉包 | rank4開のPfaffianグラフ閉包＝Rees blowup | prop:graph。純rank2成分は別 |
| 標識の有限性 | Hまたは次数/重複度付き因子標識の忘却 | prop:finite-mark, prop:factor-mark |
| rank2入口 | 飽和可積分rank2層Fと非零s∈det F、e=0,1,2 | prop:ranktwo。foliationの完全分類ではない |

新しい基本原理の証明を検討し、29項目の恒等式を独立計算した。元147成分の全証明を再検証したとはしない。

## 2. 復元から変形への橋

E=H0(Ω²(5)), B=H0(O(5)), d_Fα=F dα−dF∧α、qα=α²/2。
相対incidenceはd_Fα=0、qα∈C×F。Fの層でrank d_Fが一定ならkernel束上の二次条件になる。
rank jumpと特殊底だけにある成分を保持する。

qα=Fと規格化し、α_t=α+ta+t²b、F_t=F+tφ+t²ψとすると、
一次は α∧a=φ、d_Fa=dφ∧α−φdα。
二次は α∧b=ψ−a²/2、d_Fb=dψ∧α−ψdα+dφ∧a−φda。
この式は指定した因子変形の可解性を測り、全高次持上げの保証ではない。
次は既存TRIPLE→Ga9と別の分裂族の実係数に適用し、必要な補正が出る理由を説明する。

標識忘却は有限でも、閉像の重複・一般fiberのGalois結合・特殊底の成分は別に管理する。
x²=tは全体が既約だが幾何学的一般fiberは二点。tx=0は一般fiberだけでは見えない特殊底成分を持つ。
rank4開からのグラフ閉包を使うことと、一般的な因子の開だけへ制限することを混同しない。

## 3. 有用な障害・反例

**rank2への誤適用。** P=∂1∧(∂2+x1∂3)ではPf P=0、P#df=f divPは自動成立だがJacobi(1,2,3)=1。
modular復元にf≠0は必須である。

**因子数1。** Dredの因子数1には既約五次と5Hがある。後者は標識affine模型でZ=0と同値。
標準可換代数が零というだけで全対称性がないと結論しない。

**H+2Qの低次数primitive。** Z_Q={α:degα≤2,Qdα−2dQ∧α=0}、
P_Q(β)=Qdβ−dQ∧β、degβ≤1 と置く。既存QQ線形計算では、rank4/5の五配置とrank3頂点回避はdim19=image19、
rank3横断Q=1+u²+v²は20/19、rank3接触Q=w+u²は25/19。
無条件のexactnessは偽。今回この全行列を再実行していない。
横断のdu∧dv/QはQ=0≅Gm上の非零residue dlog aを持つ真のde Rham類。
接触の補集合はGm×A³でH²_dR=0だが、β_r=u dr−r du、r∈Sym²⟨v,z,w⟩の二次primitiveが
P_Q(β_r)=2w du∧dr−u dw∧dr−r du∧dwという低次数αを作る。余剰6は次数濾過に関係する。
全独立性・尽くしは既存の有限計算入力。これは局所的な研究課題であり全分類の前提にはしない。

**対称性と商。** 大域ベクトル場の次元は一般点での評価rankとは異なる。
任意の複素係数半単純場が代数的Gm作用へ積分するわけではない。
非線形座標変換は族の構成には使えても、PGL帰属の代用にはならない。

## 4. 原証明への主な入口

全支持の現行地図：`research/classification/KNOWN142_SEVEN_CANDIDATE_WHOLE_POISSON_COVER_AUDIT.md`。
このファイルの旧142＋7は証明の入力段階で、現在の成分数へ読み替えない。
A：HYPERPLANE_QUARTIC_COFACTOR_ENTRANCE、DOUBLE_QUADRIC_SUPPORT、RANK_THREE_QUADRIC_SUPPORT、多Hの全支持。
B：SMOOTH_QUADRIC_CUBIC_SUPPORT、DEGENERATE_QUADRIC_CUBIC_COFACTOR_REDUCTION。
C：papers/normal_actual_divisor_recognition_ja.tex、IRREDUCIBLE_NONNORMAL_QUINTIC_WHOLE_POINTWISE_COVERAGE、
NONPLANE_CONDUCTOR_POINTWISE_KNOWN_COVERAGE、PLANE_RANKZERO_CONDUCTOR_WHOLE_POINTWISE_COVERAGEとrank2版。
R2：RANK_TWO_SUPPORT、KNOWN_COMPONENT_AUDIT、独自stable_rank2の元証明。
境界：GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF、GA9_CONSTANT_COMPONENT、T37_SIX_TARGET_LEAF_OBSTRUCTIONS_SEPARATE_REVIEW。
以上はresearch/classification内（明記したpapersを除く）。今回読んだ合成ノートと、未読の全原証明を区別する。

## 5. 一次資料の確認と射程

- B. Pym, Constructions and classifications of projective Poisson varieties, LMP108(2018),573–632, DOI10.1007/s11005-017-0984-5, arXiv:1701.08852。§5.1、Propositions5.14/5.16。SNCでの対数形式・トーリック認識の射程を確認。https://doi.org/10.1007/s11005-017-0984-5
- M. Matviichuk, B. Pym, T. Schedler, A local Torelli theorem for log symplectic manifolds, arXiv:2010.08692v2。Theorem1.4をPDF本文・頁画像で確認。SNC locusと交わる40成分であり、全成分をSNCへ退化させる定理ではない。https://arxiv.org/pdf/2010.08692v2
- M. Corrêa, A. Muniz, Holomorphic foliations of degree two and arbitrary dimension, arXiv:2207.12880v5。本文Theorem4.1を確認。安定分岐に対する可積分性・uniruled性を固定Jacobian族への全点帰属と同一視しない。https://arxiv.org/html/2207.12880v5
- Stacks Project, Tag01OF (Blowing up), Rees表示と局所chart。https://stacks.math.columbia.edu/tag/01OF
- Stacks Project, Tag0553 (Irreducible components of fibres)。一般fiberと幾何学的一般fiber・基底変換を区別する参照。https://stacks.math.columbia.edu/tag/0553

必要な主張の一次本文を確認した範囲であり、全先行研究の徹底照合・新規性確認完了とはしない。

## 6. 再現

`python code/check_global_structure.py --output work/global-structure.json`。
普遍的な6係数と24一次微分を独立変数として検算し、標本値は使わない。
原稿はREADMEのLuaLaTeX二回コマンド。PDFは別フォルダに出力し、中間生成物をbundleへ入れない。
今回の生成・検証の実結果はreports/global_structure_validation.jsonを参照。
