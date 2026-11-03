# 二つの adjoint 共変量の完全な8次線形関係と一般二重写像

2026-11-03。前 Spin(5) packet の二次共変量 C,D を使い、全8次係数を計算する。**二つの4次関係と二つの新しい8次関係、一般 fiber がちょうど ±、像の一つの局所的な codimension 4 成分**を証明する。NEW25 方程式、四式の零点全体の尽くし、全 regular invariant 環の生成は主張しない。

## 四つの普遍的な関係

q₁=tr C²、q₂=tr(CD)、q₃=tr D² とする。前 packet の F₄=tr μ² に対し q₁=F₄/5、q₂=−F₄/10、q₃=F₄/10 である。従って

\[
 f_1=q_1-2q_3=0,\qquad f_2=q_2+q_3=0.
\]

今回、全16変数の全係数で次を確認した。

\[
 \boxed{f_3=\operatorname{tr}C^4+4\operatorname{tr}D^4-3q_3^2=0,}
\]

\[
 \boxed{f_4=-\operatorname{tr}C^4
       +6\operatorname{tr}(C^2D^2)
       +12\operatorname{tr}(CD^3)
       +8\operatorname{tr}D^4=0.}
\]

f₁,f₂ は β について4次、f₃,f₄ は8次である。20個の adjoint 座標についてはそれぞれ2次、4次である。C,D は一般に可換でない。特に tr(C²D²) を tr(CDCD) へ置き換えてはならない。

f₃ は前 packet の F₄²=400det C+1600det D と Newton の4×4恒等式からも導かれる。f₄ が今回求めたもう一つの primitive な8次関係である。

## 12個の8次式の関係空間を全係数で尽くす

cyclic binary trace words の順序を

\[
 CCCC,\ CCCD,\ CCDD,\ CDCD,\ CDDD,\ DDDD
\]

とする。これら六つの trace と、

\[
 q_1^2,q_1q_2,q_1q_3,q_2^2,q_2q_3,q_3^2
\]

の合計12式を用いた。C,D の全二次係数の分母を一律4で払い、全8次式を256倍した**整数疎多項式**として掛ける。6 trace の項数は 5727,7861,7105,8229,7083,5487 である。単項式は指数を9進符号化し、8次以下なので桁上がりはない。

次の七つを全単項式で直接検査した。

\[
 q_1^2-4q_3^2,\quad q_1q_2+2q_3^2,\quad
 q_1q_3-2q_3^2,\quad q_2^2-q_3^2,\quad
 q_2q_3+q_3^2,\quad f_3,\quad f_4.
\]

七つの係数関係は QQ 上独立である。他方、12式の全係数行列から非零の5×5 QQ minor を保存した。従って rank は**ちょうど5**で、全線形関係空間の次元は7、上の七つがその基底である。五つの積関係を除くと primitive な関係はちょうど二つになる。これは任意次数の非線形 ideal の尽くしではない。

既受理の dim Inv₈=5 と合わせると、この12式の span は全8次不変式空間である。今回は点への評価で関係を推測せず、整数係数の完全な同一性と QQ minor で証明した。

## 同じ一つの点での24×20線形復元

使用した QQ 点は `Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63` の sample[0] と**同一**である。新しい β 点は追加していない。以下はその点の計算を利用した周囲 W 全体の構造に関する議論で、NEW25 上の一般性を仮定しない。

一般に Ψ は5×4行列で

\[
 C=(\Psi^{\mathsf T}Q\Psi)J^{-1},\quad
 V(D)=(\Psi J^{-1}\Psi^{\mathsf T})Q.
\]

従って恒等的に

\[
 \boxed{V(D)\Psi-\Psi(J^{-1}CJ)=0.}
\]

これは行列積の結合則による普遍的な恒等式である。これら20線形式と4個の Clifford trace 線形式を、未知 Ψ の20係数について並べた24×20行列を L(C,D) とする。

同じ QQ 点で L の rank は**19**で、核は実際の Ψ の張る直線である。完全な L、非零19×19 minor、その行・列番号、残った free column、free coordinate を1にした実核を保存した。この minor の非零開集合上では、kernel line を19次 minors、または固定した19×19角の逆で**有理的に**回復できる。

さらに20個の adjoint 座標を、前 packet と同じ10個の Aᵢᵀ 基底を C,D それぞれに使って固定した。C,D の全20×136二次係数から作った20×16微分は同じ点で rank16であり、完全な非零16×16 QQ minor を保存した。単に形式的な線形復元系の解を得ただけでなく、元の写像の微分も検査している。

## 一般 fiber が ± であること

φ:W→sp₄×sp₄、β↦(C,D) とし、像の reduced closure を Z とする。W の次元は16であり、上の rank16により dim Z=16。L の19 minor が非零であることと、保存した C または D の一つの非零成分を条件として Z の非空開集合を取る。

同じ共変量対を持つ β,β′ は両方ともこの kernel line に属するから Ψ′=λΨ、従って β′=λβ である。C,D は二次で、選んだ成分が非零だから λ²=1。よって一般 fiber は正確に {β,−β} となる。微分 rank16 からこの二点は reduced である。実点の選んだ共変量成分と kernel minor を証明書に固定しているので、非空開集合の存在を抽象的な一般性に置き換えてはいない。

従って

\[
 \boxed{\mathbb C(Z)=\mathbb C(W)^{\{\pm1\}},\qquad
 [\mathbb C(W):\mathbb C(Z)]=2.}
\]

−I∈Sp₄ は W の三つの covariant index により β↦−β と作用する。φ は Sp₄-equivariant であるから、

\[
 \mathbb C(Z)^{Sp_4}=\mathbb C(W)^{Sp_4}.
\]

一般の二つの β の共変量対が同じ同時共役軌道にあれば、共役を施して対を一致させた後の違いは ± だけであり、これも Sp₄ の中心の作用である。従って一般軌道の分類は共変量対の同時共役分類へ帰着する。これは全ての特異点・退化軌道の分類、regular invariant 環の同型、写像の全域での有限性を意味しない。

## 明示的な四式の一つの局所成分

20次元 affine space sp₄×sp₄ 内に Y=V(f₁,f₂,f₃,f₄) を取る。全式は上の完全係数恒等式から Z 上で消える。同じ一つの点で、四式の **4×20 gradient の rank は4**である。非零4×4 QQ minor と gradient 全体を保存し、gradient×dφ=0 も全64係数で検査した。

従って Y はこの点の近傍で滑らかな16次元多様体である。Z は16次元既約閉集合で Y に含まれるから、この点で Z と Y の局所環は一致する。また dφ はその接空間への同型なので、β と −β の各近傍の完成写像は形式的な同型である。

以上は **Z が明示的な四式の零点の一つの generically smooth codimension 4 成分である**ことを与える。Y 全体に別の成分や非被約構造がないこと、φ の像が全域で閉であること、Z の全域の正規性、四式が Z の大域 ideal を生成することは主張しない。

## Character のみの追加値

既存の16ウェイトの整数DPと C₂ Weyl 交代和を24次まで計算し、

\[
 \dim\operatorname{Inv}_{22}=21,\qquad
 \dim\operatorname{Inv}_{24}=67
\]

を得た。全 Weyl 項を保持し、各次数の全重み数を binomial(15+d,d) と照合する。これは character 計算だけであり、22次・24次の不変式評価、基底構成、NEW25 上の rank、生成元や関係式の個数は計算していない。

## 保存と再現

* `code/q_firstjet_two_adjoint_relations_20261103_v63.py`
* `reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/{report.json,certificate.json.gz}`

証明書は12個の全整数8次多項式、7関係、rank5 minor、全復元行列、実 Ψ、全20×16微分、四式の全 gradient、各 QQ minor と24次までの character を保持する。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_firstjet_two_adjoint_relations_20261103_v63.py \
  --output work/two-adjoint-relations-replay
```

60秒・2GB上限。実行約1.71秒。新しい β 点、GB、22次以上の評価行列、NEW25 の方程式は含まない。
