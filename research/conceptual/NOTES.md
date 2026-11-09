# 大域構造研究ノート

2026-11-09 / v5。最初に RESEARCH_POLICY.md、次に GLOBAL_PLAN.md を読む。

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
| modular原始関数 | df+i_Zα=i_Jμ、f零も許す強い全Jacobi判定 | 共通TeXのprop:modular-primitive、今回独立恒等式検算 |
| 共鳴 | Z=S+Nから[S,P]=[N,P]=0、有限な単項式共鳴型 | prop:resonance。有限な固有値比とは言わない |
| 相対復元・持上げ | 固定Fのkernelに加え、因子も動く一次・二次式 | prop:relative, prop:lifting。全高次無障害は未主張 |
| グラフ閉包 | rank4開のPfaffianグラフ閉包＝Rees blowup | prop:graph。純rank2成分は別 |
| 標識の有限性 | Hまたは次数/重複度付き因子標識の忘却 | prop:finite-mark, prop:factor-mark |
| 斉次cofactor復元 | divΠ零、Π#dq=qX、β=d(qc)/5で全rank3/4/5。Xtrace零、対数対称性は5X/6 | quadric_cubic_reconstruction.tex、thm:quadric-cofactor-reconstruct |
| 全被約平面五次 | 自由Saito基底と全R6、隠れたfiber依存なし | thm:reduced-plane-whole-source、A/B/Cの三箇所で共有 |
| 磁気留数 | kerδ≅{λ:Σd_iλ_i=0}、ker/coker次元ν−1 | thm:magnetic-residue-quotient、通常の全F証明と八例の記号検算 |
| rank2入口 | 飽和可積分rank2層Fと非零s∈det F、e=0,1,2 | prop:ranktwo。foliationの完全分類ではない |

単純Hと重複Hの入口を共通TeXへ置換し主稿にも統合した。第三の二H+三次と定数modular全源を確認した。
旧29項目を再現し、新しい原始関数恒等式と全表示族を独立計算した。元147成分の全証明を再検証したとはしない。
全147 IDのactual一般元と出典はcomponent_entrances.csv。A113/B10/C18/R2六は一般元の案内である。

## 2. 復元から変形への橋

E=H0(Ω²(5)), B=H0(O(5)), d_Fα=F dα−dF∧α、qα=α²/2。
相対incidenceはd_Fα=0、qα∈C×F。Fの層でrank d_Fが一定ならkernel束上の二次条件になる。
rank jumpと特殊底だけにある成分を保持する。

qα=Fと規格化し、α_t=α+ta+t²b、F_t=F+tφ+t²ψとすると、
一次は α∧a=φ、d_Fa=dφ∧α−φdα。
二次は α∧b=ψ−a²/2、d_Fb=dψ∧α−ψdα+dφ∧a−φda。
この式は指定した因子変形の可解性を測り、全高次持上げの保証ではない。
共有TeX poisson_lifting_examples.tex にTRIPLE→Ga9と可換場の二重H分裂を適用した。
TRIPLEではdegree≤2がA=e/12,B=e(15−e)/72を強制し、q(a)=0だがψ≠0なので指定因子族に二次補正bが必要。
一次tensorだけではJacobiの(t,x,r)成分にe²(x²+5y)/12が残る。
分裂例は可換場と分解可能方向によりb=0で全次数が閉じる。任意方向への一般化はしない。

標識忘却は有限でも、閉像の重複・一般fiberのGalois結合・特殊底の成分は別に管理する。
x²=tは全体が既約だが幾何学的一般fiberは二点。tx=0は一般fiberだけでは見えない特殊底成分を持つ。
rank4開からのグラフ閉包を使うことと、一般的な因子の開だけへ制限することを混同しない。

## 3. 磁気商と留数の統合

詳細証明は共有TeXへ置く。ν個の相異なる因子を持つ被約free平面五次では、
δ(ℓ,m)=div(ℓV−mU)のkernelからη=i_E i_WΩ/Fへ写すと、次数関係を持つ対数留数空間になる。
逆はi_WΩ=d(Fη)/5で多項式として作れ、自由syzygyからℓ,mを一意に戻せる。
従ってkernelとcokernelはν−1次元。両空間そのものを無根拠に同一視しない。
ν=2で磁気一変数、非零23軌道と零22境界、ν=1で全磁気消去、ν=3では商二次元。
これはBの六配置のrank五とCのrank六の同じ通常の証明であり、単なる行列の再実行ではない。
因子が衝突し非被約になると適用できない。

全31保存源は全自由係数と全残Jacobiを保持したexact projectionで、全係数恒等式を再生成した。
その線形完全性のfull150 rank証明、非平面実変形・graph帰属、六平面配置の尽くしと周囲障害は旧入力として残る。
八例の検算は六Q+C、二H+C、既約Fの全R6・shear正逆・留数kernelを含む。全147の独立検証ではない。

## 導手上の共通復元

Cの導手支持について、rank零・二の可換余法線とflat作用を共通化した。
planeではcanonical liftの全weight源を(45,30,9)、rank零を(35,30,9)として正逆付きで復元した。
restriction・余法線作用を固定してもbinary quadratic三係数が残る。全higherと非斉次Jacobiを保持する。
leading Pfaffianをrestrictionで割ったsectionの平行輸送から、二次判別式と四次不変量を同じ原理で説明する。
rank零可換性、rank二module、既約cubic零因子の全source、m2/m4輸送の重複証明を実際に置換した。
全150係数で制約rank42、Euler24、全源84の完全性と全自由係数恒等式を記号検算した。
新規性未確認。旧有限場合分け・共鳴・monodromy・全支持帰属・成分性は保持し、全147の独立検証ではない。

共有証明はconductor_conormal_reconstruction.tex。全パラメータの抽出は定数2/3/5の除算のみで、
parameter零・共鳴の除外を伴わない。余法線作用のtraceは5T、作用kernelは三係数Kである。
作用零でもKを捨てず、全weight零が零のときだけ三十九higher源となる。
rank二開ではG_m/ρがSym^m N*⊗detNの平行sectionで、m2判別式はcρ²、c零も残す。
局所cuspの留数du/(2v)=da/a²は、接続の極から正規化後留数の単極を推測できないことを示す。
導手支持の非被約scheme構造やfamily正規化のflatnessをこの源の前提としない。

## 4. 有用な障害・反例

**rank2への誤適用。** P=∂1∧(∂2+x1∂3)ではPf P=0、P#df=f divPは自動成立だがJacobi(1,2,3)=1。
弱いmodular復元P#df=fZにはf≠0が必須である。
今回の強い原始関数式df+i_Zα=0はf零でも有効だが、この反例はi_Zαが非零なので除かれる。

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

## 5. 原証明への主な入口

全支持の現行地図：`research/classification/KNOWN142_SEVEN_CANDIDATE_WHOLE_POISSON_COVER_AUDIT.md`。
このファイルの旧142＋7は証明の入力段階で、現在の成分数へ読み替えない。
A：HYPERPLANE_QUARTIC_COFACTOR_ENTRANCE、DOUBLE_QUADRIC_SUPPORT、RANK_THREE_QUADRIC_SUPPORT、多Hの全支持。
B：SMOOTH_QUADRIC_CUBIC_SUPPORT、DEGENERATE_QUADRIC_CUBIC_COFACTOR_REDUCTION。
C：papers/normal_actual_divisor_recognition_ja.tex、IRREDUCIBLE_NONNORMAL_QUINTIC_WHOLE_POINTWISE_COVERAGE、
NONPLANE_CONDUCTOR_POINTWISE_KNOWN_COVERAGE、PLANE_RANKZERO_CONDUCTOR_WHOLE_POINTWISE_COVERAGEとrank2版。
R2：RANK_TWO_SUPPORT、KNOWN_COMPONENT_AUDIT、独自stable_rank2の元証明。
境界：GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF、GA9_CONSTANT_COMPONENT、T37_SIX_TARGET_LEAF_OBSTRUCTIONS_SEPARATE_REVIEW。
以上はresearch/classification内（明記したpapersを除く）。今回読んだ合成ノートと、未読の全原証明を区別する。

## 6. 一次資料の確認と射程

- B. Pym, Constructions and classifications of projective Poisson varieties, LMP108(2018),573–632, DOI10.1007/s11005-017-0984-5, arXiv:1701.08852。§5.1、Propositions5.14/5.16。SNCでの対数形式・トーリック認識の射程を確認。https://doi.org/10.1007/s11005-017-0984-5
- M. Matviichuk, B. Pym, T. Schedler, A local Torelli theorem for log symplectic manifolds, arXiv:2010.08692v2。Theorem1.4をPDF本文・頁画像で確認。SNC locusと交わる40成分であり、全成分をSNCへ退化させる定理ではない。https://arxiv.org/pdf/2010.08692v2
- M. Corrêa, A. Muniz, Holomorphic foliations of degree two and arbitrary dimension, arXiv:2207.12880v5。本文Theorem4.1を確認。安定分岐に対する可積分性・uniruled性を固定Jacobian族への全点帰属と同一視しない。https://arxiv.org/html/2207.12880v5
- Stacks Project, Tag01OF (Blowing up), Rees表示と局所chart。https://stacks.math.columbia.edu/tag/01OF
- Stacks Project, Tag0553 (Irreducible components of fibres)。一般fiberと幾何学的一般fiber・基底変換を区別する参照。https://stacks.math.columbia.edu/tag/0553

Gualtieri–Pym, arXiv:1203.4293v1、§4.1・Example4.4・Proposition3.11の一次PDF本文を確認した。
Poisson module、canonical moduleと退化台保存の枠組みへの参照であり、planeの全係数源を引用したものではない。https://arxiv.org/abs/1203.4293

必要な主張の一次本文を確認した範囲であり、全先行研究の徹底照合・新規性確認完了とはしない。

## 7. 再現

`python code/check_global_structure.py --output work/global-structure.json`。
普遍的な6係数と24一次微分を独立変数として検算し、標本値は使わない。
原稿はREADMEのLuaLaTeX二回コマンド。PDFは別フォルダに出力し、中間生成物をbundleへ入れない。
今回の生成・検証の実結果はreports/conductor_conormal_validation.json、Bの前段はreports/quadric_residue_validation.jsonを参照。
新検算は `python code/check_marked_reconstruction.py --output work/marked-reconstruction.json`。
組版は `bash scripts/build_structural_manuscripts.sh`。

Bと磁気留数の検算は `python code/check_quadric_residue_reconstruction.py --output work/quadric-residue-checks.json`。

Cの新検算は `python code/check_conductor_conormal_reconstruction.py --output work/conductor-conormal-checks.json`。
