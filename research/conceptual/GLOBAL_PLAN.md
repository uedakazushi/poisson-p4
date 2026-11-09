# 全分類の概念的再構成：標識付き復元と持上げ可能な変形

2026-11-09 / v69-quadric-residue-reconstruction。最初に `RESEARCH_POLICY.md` を読む。

## 1. 決定した方針

研究の主軸を、局所的なH+2Q補題から、**全点の認識 → 標識付き復元 → Poisson変形 → 極大閉像**へ移す。
目指すのは既存147項に新しい名前を付けることではなく、複数のケース別入口・計算・包含証明を共通定理で実際に置き換えることである。
完全な全軌道分類、全交叉図、全非被約局所環は要求しない。147は現行稿の結論であり、この再構成の仮定ではない。

決定的な区別は二つある。
- 固定した因子上の構造の復元と、その因子を動かすPoisson変形は別問題である。
- 幾何学的な型の分類と、族の閉包が全Poisson空間の成分であることは別問題である。

この二つをつなぐ持上げ・極大性を省いた「概念化」は採用しない。

## 2. 全点を覆う少数の入口

`papers/poisson_divisor_symmetries_ja.tex` の thm:three-entrances は次数五だけから次を証明する。

| 入口 | 全点での仮定 | 主な構造と復元 | 保持する現行の入力 |
|---|---|---|---|
| A：超平面を持つ | actual Dに少なくとも一つHがある。重複度任意、5Hを含む | Hを標識し、一意なaffine quadratic Pへ。f=Pf PとZ=div Pの復元 | H+quartic、多H+cubic、三H+Q、五H、H+2Q、重複Hの各全源 |
| B：超平面なしの可約型 | 必然的にD=Q2+C3、両因子既約・被約 | 標準対数対称性をaut(Q)内で分類し、作用に適合する閉形式を復元 | SMOOTH_QUADRIC_CUBIC_SUPPORT、DEGENERATE_QUADRIC_CUBIC_COFACTOR_REDUCTIONと全支持定理 |
| C：既約五次 | actual Dが被約既約。normal/nonnormalを分離 | normalの認識、nonnormalの正規化・導手・留数・余法線記号 | normal_actual_divisor_recognition、IRREDUCIBLE_NONNORMAL_QUINTIC_WHOLE_POINTWISE_COVERAGE |
| R2：一般階数二 | Pfaffian sectionが恒等零 | 飽和rank2 foliationとdet Fの切断。degree d+zero degree e=2 | RANK_TWO_SUPPORT、独自の安定分岐証明、因子付き分岐、Aff32のrank4変形 |

これは点の入口であって成分の分割ではない。同じ成分がA/B/Cにまたがり得る。
Dredの因子数1には既約五次だけでなく5Hもある。因子由来の可換代数が零でも、全対称性が零とは限らない。

## 3. 今回得た共通数学と、まだ得ていないこと

**証明を収録・今回検討した共通原理（新規性未確認）**

1. `thm:marked-affine`：不変Hを標識するとh-Casimir liftは一意、全affine二次Poisson模型と対応。actual D=h^5 f(x/h)、mult_H D=5−deg f。
2. `thm:modular-reconstruct`：f≠0とaffine Zを固定すると、div P=Z、P#df=fZはPについて線形、残る条件はPf P=f。これらからJacobiが従う。Z=0 iff fが非零定数 iff D=5H。
3. `prop:resonance`：Zの半単純部分に対する単項式共鳴で有限個の構成可能な型へ分けられる。ただし固有値比は一般には連続であり、translation余核と全zero complementを保持する。
4. `prop:relative`/`prop:lifting`：全因子上の相対閉形式模型と、指定した因子変形の一次・二次持上げ方程式。固定因子の線形核だけでは見えない障害を明示する。
5. `prop:graph`：rank4開のPfaffianグラフ閉包はRees blowup。rank2境界と純rank2成分を分離する。巨大なRees ideal計算は要求しない。
6. `prop:finite-mark`/`prop:factor-mark`：因子標識の有限性。幾何学的一般fiberの分裂をそのまま成分数にしない。
7. `prop:ranktwo`：rank2構造は飽和foliationと零因子から復元され、e=0,1,2で尽くされる。

今回の独立SymPy検算は29項目の普遍多項式恒等式・次数展開・rank2反例であり、全分類の検証ではない。

**残る重要な数学**：可能な(f,Z)の認識、各線形・二次fiberの全支持、全源の帰属、実変形の存在、高次障害、相互非包含、正規/非正規五次認識。復元式を別の式にしただけでこれらを解いたとは扱わない。

## 4. 実行順序：まず全体へ効く一つの置換を完成する

### P0：共通入口を主稿へ統合済み

v68では全147 IDをactual記録から一表へ対応付けた（A113/B10/C18/R2六）。
共通TeX `papers/marked_hyperplane_reconstruction.tex` を両原稿が読み込み、単純Hと重複Hの導入、
modular対称性・共鳴と最高部Koszul条件の導出を同じ証明へ置換した。
modular原始関数式はf零にも成立する強い判定であり、弱いP#df=fZとの区別を保持する。
定数modular全源は加法定数を含め復元し、第三の二H+三次模型も正逆付きで確認した。
個別source全帰属と全有限計算は保持する。以下は適用域と保持条件である。

main稿の節・命題とCOMPONENTS.jsonの全項目をローカルcloneで読み、各項目をA/B/C/R2へ対応付ける。
必要な表は一つだけ：ID、generic actual D、入口、全源定理、成分性・相違の根拠、主稿の所在。
IDからDを推測しない。不明な欄を推測で埋めない。旧1158主張の総監査を前提作業にしない。
この表は独立した研究成果ではなく、次の共通証明が何を置き換えたか確認するための道具である。

**最初の数学的成果物**は、Aのうち少なくとも「単純H+四次」と「重複H」の二つの全入口を、
一つの標識超平面・modular復元・共鳴の命題として主稿へ統合すること。
比較対象は `HYPERPLANE_QUARTIC_COFACTOR_ENTRANCE.md`、
`KNOWN142_SEVEN_CANDIDATE_WHOLE_POISSON_COVER_AUDIT.md` §4、およびそれらが指す元の全source。
第三の配置（複数H、またはH+2Q）でも同じ命題の適用条件を確認し、再利用性を試す。
全sourceの一般点だけでなく、inverse・分母・translation・Pfaffian次数低下を保持する。
具体的に不要になった旧補題・重複議論と、なお必要な有限計算を明記する。
最初のCodex作業を、この計画の言い換えや監査表だけで終了しない。

### P1：二つの実持上げ機構を主稿へ統合済み

`papers/poisson_lifting_examples.tex` はTRIPLE→Ga9の次数上限制約からA=e/12,B=e(15−e)/72を導く。
指定因子族の二次項ψは非零、一次方向は分解可能でq(a)=0なので、二次tensor補正bは必要である。
一次だけではJacobiにもe²(x²+5y)/12が残る。真の全族、全正逆、actual affine帰属を区別して保持した。
別の二重H分裂族では可換場と分解可能性から全高次補正が不要となることを証明した。
これで任意の因子変形の無障害性を示したとはしない。以下の一般fiber・閉像の課題は残る。

相対持上げ式を、既存のTRIPLE→Ga9の実曲線、少なくとも一つの別の因子分裂/衝突族に適用する。
中心の構造・一次方向・二次障害・全域の非零性がどの式に対応するかを導出する。
既存係数を再実行してPASSとするだけでなく、自由度または強制される補正項が何故現れるか説明する。
因子の変形だけからPoisson変形を推測しない。非線形変換による族の構成とPGL帰属を分ける。

復元族の一般fiberの既約成分を関数体上で調べ、標識忘却の重複・Galois軌道・特殊底成分を処理する。
全ての相互交叉は不要。全点の有限閉被覆を得た後、極大性・重複排除に必要な比較だけ行う。

### P2：BとCをそれぞれ本質的な構造定理にする

Bの共通入口・全復元をv69で統合した。`papers/quadric_cubic_reconstruction.tex` は
canonical divergence-zero liftの一意性、5β=dDによる全Jacobi、非零q-cofactorとtrace零、
canonical対数対称性との係数5/6、全Sylvester余核を証明する。
rank三・四の全31代表の自由係数・残Jacobiは元sourceからexactに投影して保持し、全係数恒等式を再生成した。
線形sourceの完全性の旧full150 rank証明と、個別の全支持帰属・零補集合はそのまま依存する。
rank五と特異qの二つの入口を共有証明へ置換した。

平面型は被約Fの全Saito源を因子分割によらず復元し、A二H+C・B Q+C・C既約五次の三箇所へ適用した。
磁気写像δ:S1²→S2のkernelは因子の次数関係を持つ留数空間に同型で、ν因子なら
kernel/cokernelはν−1次元である。余核そのものを留数空間と同一視せず、次元とkernelの同型を区別する。
二因子のrank五と既約Fのrank六の二つの旧証明を統合した。
二因子の磁気零は22次元境界、非零は23次元軌道。三因子の商次元二を捨てない。
六Q+C・二H+C・既約Fの八例で全R6、正逆、留数核の全体を検算した。
被約性を失う因子衝突へ次元式を延長しない。非平面の実変形、有限六配置の尽くし、周囲障害はなお不可欠である。

次の主作業はCの非平面導手・normal認識の接続である。
Cのnonnormalでは、導手曲面Sを一つ標識する。一般restriction rank0/2を分け、
cotangent Lie algebroidとconormal上の一次作用、留数、正規化への降下が元のcoHiggs/scalar/zero-symbol分岐をどう組織するかを具体化する。
単に「Lie algebroid」と名付けて全higher jetを捨ててはならない。
Σδ_i deg S_i=6−gという有界性の仮定を確認し、有限個の低次数曲面の認識へ進む。
主な元入力はNONPLANE_CONDUCTOR_POINTWISE_KNOWN_COVERAGE、PLANE_RANKZERO_CONDUCTOR_WHOLE_POINTWISE_COVERAGEとrank2版。
normalの楕円認識は別の完成した命題として接続する。normalization-flatnessや導手の被約性を暗黙に追加しない。

R2では既存のfoliation全支持をprop:ranktwoへ接続し、純rank2成分とrank4閉包を分ける。
Corrêa–Munizの安定分岐の可積分性だけをJacobian族への全点帰属と読み替えない。

### P3：主稿の本体を再構成する

順序は「主定理と幾何的総覧 → 全入口・標識 → 共通復元と持上げ → A/B/Cの認識 → R2 → 極大閉像と一般元表」。
S61・134式・個別minorは必要な証明を支える座標模型として保持し、読者の最初の入口にはしない。
有効な結果や十分な証明は現行稿へ統合する。旧稿を別に置くだけで保存したことにしない。
原稿を分割するかはこの依存構造を整理してから決める。既存の価値ある局所結果も無断で削除しない。

## 5. この方針で強制しないこと

H+2Qの0/1/6余剰の完全な概念的導出は局所課題であり、全体の必須前提ではない。
すべてをトーリック退化へ押し込めない。MPSの40はSNC locusと交わる成分で、全141 rank4成分ではない。
全可換Lie部分代数・全foliation・全正規化の一般論の分類へ問題を広げない。
全軌道・全scheme厚み・正確nil/H0/T2・全147²交叉表は主目標にしない。

## 6. 失敗を避ける条件と終了判定

- 固定f,Zの線形性を、f,Zも動く問題の線形性と混同しない。
- f=0ではmodular–Pfaffian判定はJacobiを保証しない（原稿に反例）。
- constant-rank kernelの記述を全底へ無条件に延長しない。縦のrank4成分を捨てない。
- rank4開での飽和と、一般因子層だけへの局所化を区別する。
- 一般点の因子型/Jordan型の相違は閉包非包含の証明ではない。
- rank3のH+2Qでaffine-linear primitiveの一様存在を復活させない。

一つの作業単位の終了条件は、通常の証明を持つ構造命題、その全適用域と例外、
実際に置換した複数の旧証明箇所、なお必要な入力が原稿に統合されることである。
二回続けて形式的な言い換えしか得られなければ、原因を数学的に診断して別の入口へ切り替える。
計算回数・小行列式・監査件数・PDF頁数の増加を進捗の代用にしない。

## 7. 現在の保存範囲

主分類稿の既存の全証明を保持し、共通入口と実持上げを共有TeXで主稿・構造稿へ統合した。
両PDFを二回組版し、今回の追加箇所を目視する。成分・主張台帳と旧証明書は未変更。
全147の独立検証又は文献上の新規性確認完了は主張しない。
現行の収録範囲・依存・保存状態はbundle manifestとvalidation reportに明記する。
