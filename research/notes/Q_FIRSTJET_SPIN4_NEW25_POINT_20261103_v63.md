# 一つの既存 NEW25 点による混合判別式候補の除外

2026-11-03。受理済み Spin(4) 混合判別式の構成を、**既存の NEW25 の QQ 点一つだけ**で厳密に検査した。結果は **D≠0** である。従ってこの混合非横断 divisor と、その有理的 Sp₄ 降下の極を払った多項式は、NEW25 hypersurface の方程式にならない。

これは候補の否定に十分な exact counterexample である。高次数不変式の評価行列、新しい点の探索、NEW25 の方程式次数への追加下界は含まない。一般の β に対する Spin(4) 二形式・混合判別式の構成自体は変わらない。

## 一つの点と元の β

`Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63/samples.json` の `actual_QQ_source_points[0]` を固定した。36点からの探索や候補点の入れ替えは行わない。この点は既存 dense NEW25 source の厳密な QQ 評価で、value pivot=2051/12、value Pfaffian=117727/2592 は既に非零と検証されている。

今回も保存された β₁₆ と元の24×16 primitive 基底から β₂₄ と Ψ₅×₄ を再構成した。最大 minors の cofactor 列 χ、normal vector n=Q⁻¹χ に対し

\[
 q=Q(n,n)=-\frac{10025632819168382287872}{1019784129324739375}
       =16\det C\ne0.
\]

従って既存候補の定義域 p₀=det C≠0 に実際に入る。

## QQ(s) の symplectic eigenframe

K=QQ[s]/(s²−q) を使う。q<0 なので QQ 上で非平方であり、K は二次体である。第四根や数値近似は使わない。Γₙ=ΣnᵢΓᵢ に対し Γₙ²=qI を全係数で照合し、

\[
 P_\pm=\tfrac12(I\pm\Gamma_n/s)
\]

の列から各2次元固有空間の二基底を選んだ。各空間の J⁻¹-pairing が非零の列対を選び、第二の列をその pairing の逆で割ることにより、基底行列 F を

\[
 F^{\mathsf T}J^{-1}F=J^{-1},\qquad
 F^{-1}\Gamma_nF=\operatorname{diag}(s,s,-s,-s)
\]

とした。使用した列番号、割る前の pairings、F、その逆転置 g=F⁻ᵀ を全て保存する。gᵀJg=J、det g=1 も厳密に検査した。

g による元の β の pullback を、線形変数の変換も含めて**全24係数**で実行した。結果は最初の Ψ 行が0で、閉形式・primitive の全8線形制約を満たす。cross 4×4 行列の determinant は s であり、normal の規格化も一致する。

変換後の12個の二形式係数は、既受理の mixed-coefficient の2倍規約に従って U(2,1),V(1,2) へ移した。

## 判別式と field involution

保存済み170-node回路を K 上で評価した。V の quadratic resultant R は非零であり、V=0 は滑らかな rational curve である。F₅=U(−v₁,v₀;b) の全6係数、判別式、R、D=Disc₅(F₅)/R² を全て厳密に保存した。

\[
 \boxed{D\ne0,\qquad D\in\mathbb Q.}
\]

完全な有理数は `certificate.json.gz` に保持する。分子は420 bit、分母は319 bitであり、補助的な還元値は D≡5237 (mod 32003)。非零の根拠は QQ 上の完全な有理計算であり、この還元値だけを根拠としていない。

全12入力に s↦−s を施して同じ回路を評価し、全170個の中間 node が元の値の共役になることを確認した。D の s 係数は厳密に0で、両埋め込みで同じ QQ 値になる。さらに Disc₅(F₅)=DR² を体の中で直接検査した。

この点では二形式の五交点はすべて単純である。p₀≠0 かつ D≠0 だから、どの k≥0 に対しても、同じ有理候補の多項式 numerator P=p₀ᵏ𝒟 はこの NEW25 点で非零になる。従って極の位数の問題を解いても、この候補が NEW25 方程式になる可能性はない。Spin(4) 構成に関する極の位数そのものや、別の判別式の同定は今回決定していない。

## 保存と再現

* `code/q_firstjet_spin4_new25_point_20261103_v63.py`
* `reports/Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63/report.json`
* 同 directory の `certificate.json.gz`：一つの source 点、全 β、normal、projectors、eigenframe、g、変換した全 β、二形式、quintic、全170 node、全 exact 出力。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_firstjet_spin4_new25_point_20261103_v63.py \
  --output work/spin4-one-new25-replay
```

60秒・2GB上限、実行約0.26秒。既存の36-point QQ 証明、Spin(5) 基底、mixed-discriminant 回路とその独立レビューの SHA256 を報告に固定した。高次数 probe や2点目の検査へ自動的に進まない。
