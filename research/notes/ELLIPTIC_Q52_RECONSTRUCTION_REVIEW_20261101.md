# q(5,2) の一次ジェットからの復元：独立点検

2026-11-01。elliptic_jets 担当の inverse-kernel 論証を独立点検した。
結論：T40=q(5,2)→W16 の一般次数は **1**。全 S61 の尽くしや T39 の次数はこの結果に含めない。

## 1. 固定された五次元核の exact check

同族のアフィン双ベクトルを P(y)、Ω(y)=−P(y)⁻¹ とし、
A_y=[−y|I₄] を V₅/⟨(1,y)⟩ への射影とする。定数 bivector φ∈Λ²V₅ に対する条件

\[
 \langle\Omega(y),A_y\phi A_y^t\rangle=0\quad\text{for all }y
\]

の核 K_P を考える。pairing は i<j の和。これは射影幾何的な定義であり、全体の
スカラーや斉次双ベクトルの radial lift の選択には依存しない。PGL₅ は K_P に自然に作用する。

`code/elliptic_q52_reconstruction_review_20261101.py` は提案された Φ_λ(r) の五つの
定数 bivector を直接代入し、Pf(P)Ω を用いて五つの恒等式を ℚ(λ)[y₁,…,y₄] で確認した。
さらに係数行列の階数5の小行列式と独立な五つの核元を確認したため、これらは全核である。
原点での spatial projection K_P→Λ²ℂ⁴ にも非零5次小行列式がある。
従って一般点では、像は Ω(0) の pairing の核である五次元超平面に等しい。

証明書：`reports/ELLIPTIC_Q52_RECONSTRUCTION_REVIEW_20261101.json`。
元の inverse-kernel 探索コードを import せず、表示された核の直接検証を行った。

## 2. 一次微分だけによる核の復元

基点を原点とし、全 Poisson 双ベクトルを既存の S61 規約で正規化する。
P(0)=J、Ω(0)=J とし、β=−dΩ|₀ を primitive にする。空間を U⊕ℂe₀ の順に並べる。
spatial projection が同型である開集合では K_P は

\[
 \phi_Z=\begin{pmatrix}Z&w(Z)\\-w(Z)^t&0\end{pmatrix},
 \qquad Z\in\ker\langle J,\cdot\rangle
\]

という一意の graph になる。点 y で射影した行列は

\[
 Z+yw^t-wy^t.
\]

\(\ell_\beta(Z)_a=\sum_{i<j}\beta_{a,ij}Z_{ij}\) と置く。核条件の y に関する
一次項は

\[
 -\ell_\beta(Z)\cdot y+y^tJw=0.
\]

従って

\[
 \boxed{w(Z)=-J\ell_\beta(Z)}.
\]

符号と係数2の有無は、pairing を i<j で固定した上で形式変数によって独立に点検した。
β は K_P 全体を一意に復元する。

補足として、閉な一次二形式20次元から graph の20係数への線形写像が同型であることも
点検した。一般24係数からの写像の核は a⊗J の四次元であり、これが閉なら a∧J=0、
すなわち a=0 となる。従ってここで失われる trace 成分が別の閉な一次微分を隠すことはない。
primitive 条件はこの20次元空間の既存の16次元 gauge section を選ぶものである。

## 3. 核から全 Poisson 双ベクトルへ

同じ primitive β を持つ、一般の正規化された二つの q(5,2) 元 P,P′ を考える。
前節により K_P=K_{P′}。一般の y で射影した K は五次元であるから、その annihilator
は一次元。従って Ω′(y)=f(y)Ω(y) となる。

両方は非退化 Poisson 双ベクトルの逆なので dΩ=dΩ′=0。よって df∧Ω=0。
四次元 symplectic 空間では α↦α∧Ω は一次形式から三次形式への同型なので df=0。
f は定数であり、Ω′(0)=Ω(0)=J により f=1。従って P′=P。

ここでは Poisson Cremona 同値が一次ジェットを保つとは仮定していない。復元は
元の射影空間にある定数 bivector の核を使う。非線形 Cremona 写像を正規化の許容
座標変換に混ぜる必要はない。

T40→W16 の優越性・一般有限性は v27 で既に証明されている。上の一般単射性と
標数0を合わせると、関数体次数は1である。全ての特殊 β 上で一意延長が存在する
という主張ではなく、明示した開条件を満たす一般点についての結論である。

## 再現

```sh
.venv/bin/python code/elliptic_q52_reconstruction_review_20261101.py \
  --output /tmp/elliptic-q52-reconstruction-review.json
```

主張の依存は、既存の q(5,2) 族同定・優越性、今回の全 λ,y に関する核恒等式、
そして上記の閉二形式の幾何学的議論である。一次資料を用いた canonical curve と
Grassmannian linear section の解釈はさらに有用だが、一般単射性の上記証明には
新しい同型類の仮定を追加していない。
