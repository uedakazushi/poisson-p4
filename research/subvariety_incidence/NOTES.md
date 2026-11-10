# Poisson部分多様体incidence：現行研究ノート

2026-11-10 / v8。証明の本文は papers/poisson_subvariety_incidence_ja.tex とその共有TeX。再開時は RESEARCH_POLICY.md を先に読む。

## 完全等方Log点での周囲二次障害

q=u0²+2u1 w0+2z w1, P0 は i_(P0)Ω=dq∧i_ξ(du0∧du1∧dz) (ξ=(1,1,-2))。
既存の固定Q・σ延長ファイバーは Spec C[a,b,c]/(ac,c²)。a=b=c=0の方向cを、C=u1∂u1−z∂z−w0∂w0+w1∂w1 として K=q(AG^-1)/2+E∧C/2 で代表する。

**新しい証明**：全ambientのPoisson条件はF(P)=E∧Jac(P)=0。一次方向KではdJac_(P0)K=0だが、二次補正Rがあるには E∧(dJac_(P0)R+Jac(K))=0 が必要。
四重ベクトルFの係数抽出
Λ(F)=2[u0³z]F0123+[u0u1zw0]F0123+[u0z²w1]F0123+4[u0³u1]F0124+10[u0u1²w0]F0124
は任意の斉次二次双ベクトルR（150係数）に対しΛ(E∧dJac_(P0)R)=0。一方Λ(E∧Jac(K))=-9。よってどのambient二次補正も不可能。5係数の縮約行列（5×6）は現行TeXに表示。全射影Jacobi線形化のrank=95、augmented rank=96、接空間次元=30は厳密なQQ計算で確認。コードと証明書：check_ambient_quadratic_obstruction.py、ambient_quadratic_obstruction.json。

**重要な論理的区別**：固定ファイバーの非被約性は証明済み。ambientの二次障害が非零なのでambient Poissonスキームはこの点で特異だが、ambientのnilradicalが非零とはまだ言えない。被約特異スキームでも接方向が二次まで延長できない例はある。

## 幾何学的な意味と残る障害

## 二次接錐の現行結果と残る問題

既存の五係数消滅汎関数はambient Poissonスキームの30次元Zariski接空間上に非零二次式を与え、完全等方Log点の局所次元を29以下に抑える。固定quadric延長の三方向(a,b,c)への制限は正確に -9c²。36次元のR2_RAT113と39次元のR2_LPBはこの点を含まない。証明はpapers/quadric_ambient_quadratic_obstruction.texの系cor:ambient-isotropic-tangent-cone-bound。三方向の直接展開は今回独立に検算した。全接錐の成分分解、局所次元の確定、固定ファイバーのacのambient障害、ambient非被約性は未解決。次はカタログの次元29以下の族と退化因子Q+2H+H'を照合し、帰属を特定する。

固定ファイバーのstrict場をA1=u1∂u0−u0∂w0、A2=z∂w0−u1∂w1、A3=Cとすると、二次式の五抽出係数は(c²/2,c²,c²,−c²/2,−c²)、重み(2,1,1,4,10)で−9c²。別実装の有理数シンボリック計算で検査した。


既存147成分分類は被約成分の分類。新しい結果はSchouten/Kuranishi obstructionが相対条件を外しても残るというscheme側の現象。次にrank4が稠密な被約A²成分が147成分のどれに属するか、P0でのtangent coneおよび各成分の接空間との関係を調べる。証明済みの主張と探索すべき仮説を区別する。

前の成果：一般Q3 Rat/Log/Aff（被約一点/A¹/被約二点）、共鳴Log非退化ブロックの直交壁、q|Wが退化する三つの滑らかなQ上の正確な全延長、滑らかなcubicの一意延長、一般del Pezzo剛性、単純HのPoisson横断対称性を現行原稿・共有TeXに保持。文献上の新規性は全面未確認。

次点は特殊Log/Aff、del Pezzo例外、P3の六成分の横断対称性。有限場合分けが価値を下げれば構造定理やより有望な応用へ切り替える。

GitHub mainだけが正本。9本のproducer/9件の証明書、TeX、再生成手順を揃え、PDF非追跡、Actions不使用、Libraryは実push失敗時のみ緊急退避。
