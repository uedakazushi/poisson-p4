# 全member退化CI22 pencil：独立網羅性・全源確認

2026-11-07。producer は root、reviewer は plane_collision_review。
`ALL_DEGENERATE_CI22_ABELIAN_CONORMAL_EXCLUSION.md` を独立に確認した。
**判定：integral・非退化CI22の全member退化pencilに対する、rank零・既約被約quintic枝の排除としてPASS。**
既受理の regular-pencil V2 と合わせ、全 integral 非退化 CI22 rank零枝の G 像はdimension≤24。
CI22上の全Poisson tensorが零、可約D源の全E帰属、又はG全体の完了は結論しない。

## 1. canonical theoremは任意のpencilを扱う

[De Terán–Dmytryshyn–Dopico, *Generic symmetric matrix pencils with bounded rank*](https://gauss.uc3m.es/fdopico/papers/jst2018.pdf)
の著者PDFを直接開き、§2、printed pp.4–5のblock定義とTheorem2.1を読んだ。
同定理の対象は**every symmetric pencil**であり、complex congruenceによる
M_d、有限固有値及び無限固有値のregular symmetric Jordan blockの直和である。
M_dのsizeは2d+1、M0は1×1零block。
論文題名にあるgeneric bounded-rank closure theoremを全pencilへ転用したわけではない。
これは引用されたThompsonのcanonical theoremを外部入力として用いる箇所であり、
その原典証明全体を再証明したとは主張しない。

size五で少なくとも一つのsingular blockがある全分割を列挙すると九つとなる。
M2単独、M1+regular2、M0+regular4、M0+M1+regular1、
及びM0を二つ以上持つ場合だけで尽くされる。
regular部分の固有値を∞から外すpencil基底の変更も許すので、∞固有値を落としていない。

## 2. 有限の全型排除とcommon radicalの残余

integral非退化Sを含む非零quadricはrank≥3。
rank≤2なら複素数体上でlinear factorを持ち、素ideal I_Sにその一つが入ってしまう。

- M2標準形はplaneを含む。二次元integral CIが同次元planeを含めばそのplane自体となり、矛盾。
- M1+regular2の相異二固有値は、pencil PGL2とx1,x2の同時変更で
  `q1=x0x1+x3²,q2=x0x2+x4²` となる。
  同一固有値のdiagonalizable型はrank2 memberを持ち、Jordan2型はplaneを含む。
- M0+M1+regular1は、regular scalar部分を消すmemberが非零rank2となる。
- common radicalが二次元以上なら残るmatrix sizeは三以下。
  size三のhomogeneous determinantは零又はP1上にrootを持ち、非零rank≤2 memberを生む。

M0+regular4では一般member Q を非退化に取れる。
rank条件からself-adjoint Mの各固有値のeigenspaceは一次元であり、各primary partは単一Jordan block。
二固有値以上なら `[T,M]∈span(I,M)` の各primary block traceで右辺係数が零になり、
異なる固有値間のintertwinerも零である。各一blockのcommutantはself-adjoint polynomialだけなので、
Q-skew centralizerは零。従って全base bivectorは零で、vertex項のみのrank≤2になる。

一固有値なら単一Jordan4で、skew normalizerはgrading一次元のみ。
全sourceは `Π=alpha*A+∂c∧V`、alphaはq1,q2の一次結合である。
外積平方はalphaで割れるので、実Dもquadric alphaを因子に持つか零。
ここではregular-pencilの全5次元定理を直接使うのではなく、同じcommutator証明を四次元で適用した。
独立checkerはその全16×8 normalizer matrixのrank7・完全kernel一次元も再構成した。

## 3. M1+相異二固有値モデルの全二十源

このモデルは実際にintegral CIである。
x0≠0ではx1,x2を消去でき、irreducibleでsmoothな稠密openとなる。
x0=0のboundaryはlineなので別のsurface componentは無い。
CIのunmixed性とgeneric reduced性によりidealはprimeとなる。

Koszul列をE=Λ²TP4とtensorして、H0(E(-4))=H1(E(-4))=0、
H0(E(-2))=Λ²C5を使えば、全sourceは一意に `q1 A+q2 B`、二十係数となる。
これらの消滅は通常・外積Euler列とO(-1),...,O(-4)のcohomologyで確認でき、Sのsmoothnessを使わない。
checkerは実150係数でsource rank20、Euler-vertical rank24、合計rank44も確認した。

既受理のmultiplicity二・Frobenius論証から、generic conormalはabelianである。
prime CIではgr_I R=(R/I)[t1,t2]がdomainなのでI²はI-primary。
従ってgeneric symbolic-square条件を全ordinary I²へ移せる。
I/I²の自由性から、A,Bの各々について `dq1 A dq2∈span(q1,q2)` が必要となる。

全十skew係数とspan係数二つの15×12行列を独立に作り、rank10、完全kernel二次元を得た。
残るskew entryは正確にA13,A24で、span係数は零。
A,B双方を合わせた30×24行列はrank20であり、全必要源は

`Π=(a q1+c q2)∂1∧∂3+(b q1+d q2)∂2∧∂4`

の四変数で尽くされる。
全actual anticanonical係数をFractionで直接計算し、

`D=−x0*(a q1+c q2)*(b q1+d q2)`

との恒等一致を確認した。Jacobiを課す前から既約被約quinticではない。
因子が零となる全境界もD=0として保持する。

## 4. 正確な到達点と保存物

全member退化pencilに残る未処理型は無い。
別作者の `REGULAR_PENCIL_CI22_ABELIAN_CONORMAL_REVIEW.md` と
`CLASSIFICATION_REGULAR_PENCIL_CI22_REVIEW_V2_20261107.json` を既受理入力として接続すると、
全integral非退化CI22 rank零枝のG非優越が得られる。
regular-pencilの旧未受理probeを、この接続へ代用していない。

- `code/classification_all_degenerate_ci22_review.py`
- `research/classification/jobs/all_degenerate_ci22_review.json`
- `reports/CLASSIFICATION_ALL_DEGENERATE_CI22_REVIEW_20261107.json`

標準ライブラリFractionと既受理reviewerの疎多項式primitivesのみを用い、producer/CAS import無し。
事前gate後0.012秒で終了。root証明の読解時snapshot、全有限行列、kernel、全実装を報告に保存した。
