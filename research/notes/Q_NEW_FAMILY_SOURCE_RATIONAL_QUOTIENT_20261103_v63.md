# NEW25 source kernel の五つの第一積分と有理五変数回路

2026-11-04。前 packet の実際の30×26 raw-jet kernel を積分し、**source の微分定数体を正確に五変数の有理関数体として決定**した。明示的な rational fiber parameter と両向きの逆変換を保持し、さらに五自由変数の symbolic section から primitive β、N₁₂,m を出力する1032-node回路を保存した。

これは source の一意な冗長方向の quotient である。**actual NEW25 の Sp₄ 軌道不変量体との有限指数は未決**であり、両体が等しい、actual quotient が有理である、NEW25 方程式を得た、とは主張しない。新しい β 点・高次数 fitting は一切使わない。

## 実際の kernel の座標変換

前稿で全30 polynomial identity を検証した kernel source derivation を 𝒱 とする。source 変数は u,v,c,X,W,B、scalar 成分の表記は c∂c である。v≠0の一般部で

\[
 \delta=(v/3072)\mathcal V
\]

と rescale する。同じ一次元 foliation を与える非零有理倍である。A=256 とし、

\[
 \xi=384X,\quad\zeta=3072W-1,
 \qquad r=v(\xi+\zeta)+A,\quad s=v(\xi-\zeta)+A.
\]

ここで r,s は source-flow の座標であり、spectral eigenvalue を意味しない。逆式は

\[
 X=(r+s-2A)/(768v),\qquad
 W=\bigl(1+(r-s)/(2v)\bigr)/3072.
\]

保存済みの polynomial kernel を直接代入し、全係数で次を検査した。

\[
 \delta u=0,\quad \delta r=r(r^2-A^2),\quad
 \delta s=s(s^2-A^2),
\]
\[
 \delta v=v(r^2+s^2)+16u(r-s)^2,
\]
\[
 \delta c/c=2A^2-r^2-s^2-32u(r-s)^2/v,
\]
\[
 \delta B/B=2A^2-r^2-s^2-48u(r-s)^2/v.
\]

前稿の rank と異なる vector field を仮定していない。derivation の rescaling、全6成分、逆座標変換を記号的に照合する。

## 五つの有理第一積分

\[
 z_r=1-A^2/r^2,\quad z_s=1-A^2/s^2,
 \qquad V=v-8u(r-s)^2/(rs).
\]

上式から

\[
 \delta z_r=2A^2z_r,\quad \delta z_s=2A^2z_s,
 \quad\delta V=(r^2+s^2)V.
\]

従って

\[
 \boxed{U=u,\quad K=z_r/z_s,\quad
 I_v=V/(rsz_r),\quad
 I_c=cv^2/(rsz_r^2),\quad
 I_B=Bv^3/((rs)^2z_r^3)}
\]

は第一積分である。たとえば δlog(rs)=r²+s²−2A²、δlog(cv²)=r²+s²+2A²、δlog(Bv³)=2(r²+s²)+2A² なので、最後の二つの重みも一致する。実装は五つの δI=0を全有理式で直接検査する。

この議論の開集合では、v,r,s,zᵣ,zₛなどの表示された分母を可逆にする。元のsourceの一般部はこれらの零点に含まれない。境界を含む全域で regular な商を構成したとは言わない。

## 有理逆変換と定数体の完全性

単に五つの独立な first integrals を得ただけでは、定数体の全体を決めたことにならない。ここでは逆変換も与える。

K=zᵣ/zₛ は

\[
 (1-K)(r/A)^2+K(r/s)^2=1
\]

という conic の式である。この conic は点 (r/A,r/s)=(1,1) を持つ。line parameter

\[
 t=\frac{r/s-1}{r/A-1},\qquad D=1-K+Kt^2
\]

を用いて

\[
 R_0=1-2(1-K+Kt)/D,\quad
 S_0=1-2t(1-K+Kt)/D,
\]
\[
 r=AR_0,\qquad s=AR_0/S_0
\]

と有理的に parametrise できる。残りは

\[
 v=I_vrsz_r+8U(r-s)^2/(rs),
\]
\[
 c=I_crsz_r^2/v^2,\qquad
 B=I_B(rs)^2z_r^3/v^3
\]

と、先の X,W の式で戻る。conic、t の復元、r,s の逆合成、v,c,B の逆合成を全て厳密に検査した。従って

\[
 \boxed{\mathbb Q(u,v,c,X,W,B)
        =\mathbb Q(U,K,I_v,I_c,I_B)(t).}
\]

これは六つの新座標の algebraic independence も同時に証明する。また

\[
 \delta t=\frac{Ar(r-s)(s-A)}{s(r-A)}\ne0.
\]

従って δ は右辺で非零有理係数×∂t である。標数0で一変数有理関数の微分が零なら係数体の元だから、

\[
 \boxed{\ker\delta=\mathbb Q(U,K,I_v,I_c,I_B).}
\]

ここで rational なのは、source の冗長方向の一般 fiber と定数体である。先に構成した universal spectral-base fiber の種数を決定する結果ではない。

## 簡単な symbolic section と五自由変数 map

t=2 は、元の座標で驚くほど簡単な条件になる。

\[
 rs=A(r+s)/2
 \quad\Longleftrightarrow\quad
 v\bigl((384X)^2-(3072W-1)^2\bigr)+256\cdot384X=0.
\]

従って

\[
 \boxed{v=-\frac{98304X}{147456X^2-(3072W-1)^2}.}
\]

u,c,X,W,B の五つを自由に残し、この v を前稿の六変数回路へ合成した。t=2 は **有理な fiber parameter を固定した全パラメータ族**であり、新しい数値 β 点を評価したものではない。

全1025-node回路の各 input を明示的に置き換え、1032-nodeの `source-five-circuit.json.gz` を得た。出力は全β₁₆、代表 C,D、F₄、全Inv₈基底、p₀,…,p₄、N₁₂,m。scalar invariant 出力は source kernel の定数なので、t を失っても同じ invariant 関数を与える。βや C,D の行列そのものが t に依存しないとは主張しない。

この回路の domain は元の uvcB(32u+v) と新しい表示分母が非零の開で、m には追加で δ₈≠0を要する。全て非空な一般開集合であり、特異境界の分類をここへ含めない。

## actual quotient への有限な重複と既知 involution

前稿の raw-jet kernel は source 微分と全 group nuisance が打ち消し合う証明書である。primitive 正規化後の全 Sp₄ rational invariant は source derivation δ に殺される。従って actual NEW25 の rational orbit-invariant field L について

\[
 L\subset\mathbb Q(U,K,I_v,I_c,I_B)
\]

を得る。actual quotient は既証5次元なので、この包含は有限な代数拡大を与える。しかしその指数は今回決めていない。五変数 circuit は actual NEW25 orbit quotient への **generically finite dominant source**であり、birational 同定ではない。

少なくとも一つの有限な重複は明示できる。元の有限 parameter 変換 ε=−1,a=−1 から

\[
 (u,v,c,X,W,B)\longmapsto(-u,-v,-c,-X,W,B)
\]

を得る。コードは G=diag(−1,1,−1,1) による全 raw value と全 linear jet の有限変換恒等式として検査した。従って全 actual Sp₄ invariant を保存する。

この involution は r,s を交換し、

\[
 (U,K,I_v,I_c,I_B)\longmapsto
 (-U,1/K,-KI_v,-K^2I_c,-K^3I_B),
\]

t↦t/(t−1) と作用する。t=2 section は保たれる。section の五自由座標では (u,c,X) の同時符号反転で、W,B は固定される。よって既知 involution の固定体は

\[
 \mathbb Q(X^2,u/X,c/X,W,B),
\]

という有理関数体である。X をその平方根として戻せば二重拡大なので、固定体の尽くしも分かる。

従って L は少なくともこの既知固定体に含まれ、full source constant field からの指数は少なくとも2である。他に finite identification がないことは未証明であり、既知固定体=Lとも主張しない。

## 新しい有限問題

N₁₂ と spectral base の invariant を、actual 族から五自由変数の有理回路へ結び付けた。これにより今後は6 source invariant 間の relation を、この明示的な五変数 source、またはその既知 involution quotient から調べられる。

N₁₂に二次の関係があること、最小次数24の NEW25 方程式が存在することは仮定しない。ambient residual curve の m-degree≤8は、NEW25 multisection の base 上の次数を上から抑えるものではない。

## 保存と再現

* `code/q_new_family_source_rational_quotient_20261103_v63.py`
* `reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/{report.json,certificate.json,source-five-circuit.json.gz}`

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_source_rational_quotient_20261103_v63.py \
  --output work/new25-source-five-replay
```

上限60秒・2GB、実行約0.81秒。全 derivation、五積分、両向きの inverse、非零 δt、symbolic section、有限 raw-jet involution と合成回路を保存する。source の定数体の完全性と、actual orbit field の未知有限指数を区別する。
