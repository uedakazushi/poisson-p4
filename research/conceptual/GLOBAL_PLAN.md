# 全分類の概念的再構成：標識付き復元と持上げ可能な変形

2026-11-09 / v67-global-reconstruction。最初に `RESEARCH_POLICY.md` を読む。

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

### P0：短い全体対応付けと、一つの実装済み置換

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

### P1：復元から成分へ進む共通の変形機構

相対持上げ式を、既存のTRIPLE→Ga9の実曲線、少なくとも一つの別の因子分裂/衝突族に適用する。
中心の構造・一次方向・二次障害・全域の非零性がどの式に対応するかを導出する。
既存係数を再実行してPASSとするだけでなく、自由度または強制される補正項が何故現れるか説明する。
因子の変形だけからPoisson変形を推測しない。非線形変換による族の構成とPGL帰属を分ける。

復元族の一般fiberの既約成分を関数体上で調べ、標識忘却の重複・Galois軌道・特殊底成分を処理する。
全ての相互交叉は不要。全点の有限閉被覆を得た後、極大性・重複排除に必要な比較だけ行う。

### P2：BとCをそれぞれ本質的な構造定理にする

Bでは非零対称性X=π#((1/2)dlog q−(1/3)dlog c)を入口にし、aut(Q)の表現と閉形式復元を接続する。
Aの共鳴の考え方を再利用するが、不変Hがないのにaffine二次模型へ移さない。

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

今回の構造原稿は旧7頁原稿の有効な全証明を保持して拡張した。368頁主稿と147成分台帳は未変更。
新原理の証明・記号恒等式を確認したが、全147の源の認識・成分性・相違を独立再検証したわけではない。
完全な原稿・証明資産はGitHub上にある。この環境では1.25MBの主稿と全証明書をコンテナへ取得できず、
今回の配布bundleは**大域構造研究を再開する一式**であって、全147成分の証明bundleを代替しない。
Codexのclone環境では現行主稿と必要な依存閉包を含む完全bundleを作成し、旧ZIPや履歴は入れない。
