# Poisson部分多様体のincidence：現行研究ノート

2026-11-10 / v1

公表原稿：`papers/poisson_subvariety_incidence_ja.tex`。以下は未解決の問題と後で役立つ区別を整理したもので、論文の証明を重複させない。

## 1. 延長incidenceと幾何学的な分岐

固定したX⊂P4に対し、Poisson ideal条件で定まるambient incidenceの非零制限部分を考える。rank4 ambientでXが真のPoisson部分多様体ならXは実Pfaffian五次Dに含まれる。rank2 Jacobian型ではPfaffianが恒等零で、この制約はない。

滑らかなcubic Xでは制限写像の全fiberが一点となるため、rank4でsmooth cubic factorを持つ枝は空。quadricでは同じintrinsic bracketにrank2/rank4が共存し得る。CI22では一般のpairでrank2の一点のみであり、特殊locusが核心となる。

## 2. 一般del Pezzo剛性の証明書と障害

全ambient extensionはJ_h+q0 A+q1 B、A,B∈Λ²C5、計20変数。正規化q0=Σx_i²、q1=Σλ_i x_i²、λ=(0,1,2,u,v)の稠密chartで150係数式を作る。二次部分は符号character blockに分かれ、Q(u,v)上のrank60、90個の普遍syzygyを持つ。線形部分は(u,v,h)=(3,4,1,…,1)でrank18、残二変数の線形・二次単項式5個の係数rank5を得る。二つのnonzero minorから非空開で理想がmaximalとなる。これは全smoothを証明しない。特殊h、λのcollision、余法線作用のjumpを独立に調べる必要がある。

## 3. 次元の忘却

cubic：P34×Gr(2,5)の40次元からF3↦F3+Sym³⟨L,M⟩の4次元を引き36。
CI22：Gr(2,15)の26次元＋P4の4次元からq_i↦q_i+a_i h²の2次元を引き28。
ファイバーの開条件と一般中心の回収は論文を参照。

## 4. rank4の単純超平面

横断的なPoisson対称性vは [v,σ]=0、v∧σ≠0。これが単純H延長の必要十分条件であり、重複Hの分類と同一視しない。特別なquadric pencilとhの位置を変えると、この対称性のjumpと延長rankのjumpが関係する可能性がある（未証明）。

## 5. 既存資料と新規性

cubic/quadric intrinsic分類：LPT（arXiv:1112.3871）、Mayanskiy（arXiv:1303.5916）。smooth quadricの正規Hamiltonian全源とCI22 rank0/rank2の排除は既存のpoisson-p4研究内にある。今回の一般CI22延長剛性の文献上の優先権は未確認。計算は147成分の全再検証ではない。
