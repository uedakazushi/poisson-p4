# primitive first-jet の Spin(5) 模型と二つの adjoint 共変量

2026-11-03。新しい点の評価や22次以上の不変式探索を行わず、受理済みの16次元表現を Clifford 模型へ同定する。二つの二次共変量から作る4×4行列 pencil の行列式係数は、**8次不変式全体の基底**になる。これは周囲の表現に関する恒等式であり、NEW25 hypersurface の方程式・判別式候補の同定ではない。

## 1. 元の表現と Clifford 積

`weight1_structure.py` の規約を使う。4次元空間の交代行列を

\[
 J=\begin{pmatrix}0&1&0&0\\-1&0&0&0\\0&0&0&1\\0&0&-1&0\end{pmatrix}
\]

とする。閉じた primitive な線形2形式 β の24係数に対する8個の独立な線形制約の核を W とする。従って dim W=16。受理済みの基底 B は、今回も元の制約と全10個の物理的 Lie 微分に照合した。

primitive 定数2形式の基底 α を、係数の順序

\[
 (01),(02),(03),(12),(13)
\]

で取る。最初の元だけ α₀₁=1, α₂₃=−1 とし、他の元は表示された交代係数を1とする。Γᵢ=αᵢJ⁻¹ とすると、全25組について

\[
 \Gamma_i\Gamma_j+\Gamma_j\Gamma_i=2Q_{ij}I_4,
 \qquad Q_{ij}=\tfrac14\operatorname{tr}(\Gamma_i\Gamma_j),
\]

\[
 Q=\begin{pmatrix}
 1&0&0&0&0\\
 0&0&0&0&1/2\\
 0&0&0&-1/2&0\\
 0&0&-1/2&0&0\\
 0&1/2&0&0&0
 \end{pmatrix},\qquad \det Q=1/16.
\]

β を5×4係数行列 Ψ として表す。全16基底について

\[
 \sum_i\Gamma_i\Psi_{i,*}^{\mathsf T}=0
\]

を検査した。20次元空間から4次元空間へのこの Clifford 積は全射であり、B の像の次元は16だから、**W は gamma-traceless vector-spinor の全空間**である。ここで Spin(5)≅Sp(4) の名称はこの明示的模型を指す。

より正確には、A∈sp₄ に対し R₅(A) を α↦Aᵀα+αA の5次元行列とする。全10生成元と全16基底で

\[
 R_5(A)^{\mathsf T}Q+QR_5(A)=0,\qquad
 [A^{\mathsf T},\Gamma_i]=\sum_jR_5(A)_{ji}\Gamma_j,
\]

\[
 \Psi(\rho(A)b)=R_5(A)\Psi(b)+\Psi(b)A
\]

を検査した。元の pullback 行列 ρ は Lie 反作用であり、spinor 側は Aᵀ である。この転置を省くと moment 行列との比較の符号が変わる。

## 2. 二つの二次 adjoint 共変量

次を定義する。

\[
 C=(\Psi^{\mathsf T}Q\Psi)J^{-1},\qquad
 V=(\Psi J^{-1}\Psi^{\mathsf T})Q.
\]

C は交代形式 J⁻¹ の sp₄ に、V は so(Q) に属する。R₅ の10×10非零 minor の逆を使い、V=ΣdᵢR₅(Aᵢ) と**全25係数**で書き戻し、

\[
 D=\sum_i d_iA_i^{\mathsf T}
\]

とする。C,D は同じ4次元空間上で共役変換する二次共変量である。

全136個の二次単項式に対する4×4係数を保存した。二つの係数列には非零2×2 minor があり、独立である。さらに C₂ Weyl 交代和から

\[
 \operatorname{Sym}^2W
 =V_{(4,2)}\oplus V_{(3,1)}\oplus 2V_{(2,0)},
 \qquad 136=81+35+2\cdot10
\]

を得る。添字は ε 座標であり、W の最高ウェイトは (2,1)=ω₁+ω₂ である。従って C,D は二次 adjoint 共変量の全2次元空間の基底になる。

## 3. 既存 moment・θ との厳密な比較

既存の16次元交代形式 Ω と trace Gram を使う規約

\[
 \mu_i(b)=\tfrac12b^{\mathsf T}\Omega\rho(A_i)b,
 \qquad \mu=\sum_i(\mathrm{Gram}^{-1}\mu_*)_iA_i
\]

に対して、全16×136係数の恒等式は

\[
 \boxed{\mu^{\mathsf T}=-C+2D}.
\]

また全100組で tr(ρ(Aᵢ)ρ(Aⱼ))=12tr(AᵢAⱼ) である。従って既存の θ=ρ(μ) に対して tr θ²=12tr μ²。この比較は既存の θ/H trace 不変式と今回の模型を同じ原表現へ結び付ける。16×16の H 自体を4×4の D と同一視してはいない。

F₄=tr μ² と置く。全四次係数辞書を実際に掛けて比較し、

\[
 \operatorname{tr}C^2=F_4/5,\qquad
 \operatorname{tr}(CD)=-F_4/10,\qquad
 \operatorname{tr}D^2=F_4/10
\]

を得た。従来の `unique_quartic_sparse` の規約では F₄=(5/2)I₄ である。

## 4. adjoint pencil と全8次不変式

\[
 P(t,u)=tC+uD,\qquad
 \det P(t,u)=\sum_{i=0}^4p_i t^{4-i}u^i
\]

とする。各 pᵢ は8次 Sp₄ 不変式である。4×4行列式を置換和で**全単項式について**計算した。p₀,…,p₄ の項数は順に 2992,7323,8326,6925,5315。係数行列の5×5 minor は −1/16777216 であり、五つは QQ 上独立である。Weyl character による dim Inv₈=5 と合わせ、**五つの pᵢ は全8次不変式の基底**となる。

全係数で確認した恒等式は

\[
 \boxed{F_4^2=400p_0+1600p_4}.
\]

sp₄ 行列の特性多項式が偶数であることと上の trace 恒等式から

\[
 \det(\lambda I-P(t,u))
 =\lambda^4-\frac{F_4}{20}(2t^2-2tu+u^2)\lambda^2
   +\sum_{i=0}^4p_i t^{4-i}u^i.
\]

これで4次元 adjoint pencil と binary quartic が、座標選択を明示した普遍的な構成になった。ただしその binary quartic 判別式や、特性多項式の判別式の因子を NEW25 の方程式とする根拠はまだない。今回は判別式の展開・因数分解・NEW25 への標本評価を行っていない。

## 5. 不変式環から分かる適用限界

16ウェイトによる整数DPと全8項の C₂ Weyl 交代和を20次まで再計算した。非零の次元は

|次数|4|8|12|14|16|18|20|
|---|---:|---:|---:|---:|---:|---:|---:|
|Inv の次元|1|5|9|3|20|9|36|
|最小斉次生成元の個数|1|4|4|3|1|6|0|

16,18,20次の低次数積の rank は、既に保存・受理済みの評価行列の該当列から再抽出した19,3,36である。**新しい点への評価ではない。** Ω・表現行列・Gram 逆の全係数は32003で還元できる。従ってこれらの modular minor は QQ 上の独立性を保証し、上の character 上界と一致する。得られた生成元数の自由可換単項式数は20次まで Hilbert 次元と一致するので、この範囲には生成元間の関係式がない。全次数で自由であるとは述べない。

既存の全周囲 W の orbit minor は一般 orbit 次元10を証明している。この既存 minor の非零性を再確認し、不変式商の次元は高々6と分かる。一方、12次までに既に9個の最小生成元が必要である。よって不変式環は多項式環ではない。有限反射群の Cartan 制限によって多項式不変式環を得る通常の polar/Vinberg 模型を、そのままこの表現へ適用する道は使えない。これは新しい分類定理の主張ではなく、提案する模型が満たすべき代数的な制約である。

二つの adjoint と binary quartic は全8次を整理するが、12,14,18次に別の最小生成元が必要である。従って pencil の係数だけで全商や NEW25 方程式が決まるとも結論しない。既存の NEW25 方程式次数下界22はそのままであり、方程式自体は未解決である。

## 再現と保存

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_firstjet_spin5_structure_20261103_v63.py \
  --output work/firstjet-spin5-replay
```

60秒・2GB上限。`certificate.json.gz` は全 Clifford 行列、全10次元作用、kernel embedding、二つの全二次係数、moment の全係数、全四次 trace 辞書、pencil 行列式の全8次係数、独立 minor、character、既存積 rank の選択行列を保存する。8次単項式の疎符号は `Σ exponent_i*9^i` で、全指数は8以下だから桁上がりはない。報告は使用した元表現・既存評価 packet・コードの SHA256 を保持する。

今回の恒等式は全て QQ 上の係数恒等式である。低次数積の rank のみ既存 modular 証明を再利用する。NEW25 の式、新しい成分、判別式の因子同定、高次数探索は主張しない。
