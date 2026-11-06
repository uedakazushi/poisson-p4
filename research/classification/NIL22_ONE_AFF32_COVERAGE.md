# nil22 A=1 全源の Aff(3,2) 成分への被覆

`MODULAR_NIL22_TOP_SUPPORT.md` の A=1 全源は、既存の23次元成分
T19 に入る。以下は真の多項式族の極限による証明であり、有限 jet の積分や
局所唯一成分性を仮定しない。新成分は追加しない。

## 全源と適用範囲

`CLASSIFICATION_NIL22_ONE_CANONICAL_20261107.json` は元の七変数源全体を
次の二変数正常形へ運ぶ。座標を `(h,x,y,z,w)` とし、h は Casimir とする。

\[
\begin{array}{lll}
 P^{xy}=2hx-y^2,&
 P^{xz}=a(hx-y^2/2)+2xw+e hy,&P^{xw}=yw,\\
 P^{yz}=e h^2+yw,&P^{yw}=hw,&P^{zw}=w^2+a hw.
\end{array}                                                     \tag{1}
\]

階数四条件は `a≠0`。e は任意であり、e=0 も以下のパラメータ写像に直接含まれる。
この全源への正規化と、元の固定代表 `(a,e)=(1,159)` への24倍の関係は
上記証明書に依存する。今回の極限証明書でもその固定代表との全150係数一致を再確認した。

## 既存の実族

v14 の全族を
\[
 \Pi_{b,t}=N\wedge H_b+t\eta_b+t^2\theta_b
\]
と書く。旧座標 `(h,y,x,w,z)` では
\[
 N=h\partial_y+y\partial_x+w\partial_z,
 \quad H_b=y\partial_y+2x\partial_x+bw\partial_w+(b+1)z\partial_z.
\]
η,θ の全係数は `sources/v14.zip` の
`sections/rank2_global_eight.tex` と
`code/rank2_lower_bounds/verify_aff32_deformation.py` にある。
今回の証明書にも全十括弧を保存し、普遍的な射影 Jacobi 恒等式を再確認した。

この実族の PGL 飽和閉包が23次元全成分であり、SNC 点を含むという既存結果を使う。
正確には旧証明書は `(b,t)=(7,0)` における全 D の階数102、実族と鉛直像の
階数48を与える。これにより既約な飽和閉包は射影次元23の全成分となる。
同じ全族の `(7,1)` に SNC への一パラメータ PGL 極限があり、MPS 成分へ同定される。
台帳ラベルは T19。この成分証明と SNC 証明は今回再実行していない既存依存であり、
`audits/2026-10-30/v14_recheck/verification_summary.json` の
`aff32_universal_rank_four_deformation` に旧厳密再現記録がある。
単に rank2 Aff(3,2) 境界を共有するという理由で全族所属を推論するものではない。

## 全係数の極限

ε を新しいパラメータとして
\[
 b=-1+\varepsilon B,\qquad
 (h,y,x,w,z)_{\rm old}=(H,\varepsilon Y,\varepsilon^2X,W,H+\varepsilon Z)
                                                               \tag{2}
\]
と置く。t は固定する。この線形写像の行列式は `ε⁴`。
ε≠0 で座標変換した Π を ε 倍した族を取る。

全150係数の分母は正確に消え、ε の多項式となる。
その定数項は次の手展開でも読める。
\[
 \varepsilon N\wedge H_b\longrightarrow N_0\wedge H_0,
\]
\[
 N_0=Y\partial_X+H\partial_Y+W\partial_Z,
 \quad H_0=2X\partial_X+Y\partial_Y+BH\partial_Z-W\partial_W,
\]
\[
 \varepsilon\eta_b\longrightarrow
 (Y^2-2HX)\partial_X\wedge\partial_Z-3HW\partial_Z\wedge\partial_W,
 \qquad \varepsilon\theta_b\longrightarrow0.                    \tag{3}
\]
ここで各項は (2) による tensor 変換後の意味である。
極限の符号を反転すると六係数は
\[
\begin{array}{lll}
 2HX-Y^2,&2XW-BHY+t(2HX-Y^2),&YW,\\
 YW-BH^2,&HW,&W^2+3tHW.
\end{array}
\]
最後に `Z'=Z+tY` と置けば (1) に完全一致し、パラメータ写像は
\[
                 a=4t,\qquad e=-B.                             \tag{4}
\]
従って (1) の全 a≠0,e を直接得る。退化五次式は
`2t H²W(2HX−Y²)`。因子型は補助確認であり、所属証明には (2)–(4) の
全150係数恒等式を使う。

ε≠0 の各点は閉成分 T19 内である。全係数が ε=0 へ延長するので、
その極限も T19 に入る。(4) の全射性から A=1 の全源を被覆する。

## 証明書と独立確認

- 実装：`code/classification_nil22_one_curve.py --aff32-limit`。
- 有限契約：`research/classification/jobs/nil22_one_aff32_degeneration.json`。
- 証明書：`reports/CLASSIFICATION_NIL22_ONE_AFF32_DEGENERATION_20261107.json`。
- 今回の厳密確認：全旧族 Jacobi、全150係数の ε 多項式性、全極限、(4) の
  普遍恒等式、極限の全350 Jacobi 式、Pfaffian、元の固定代表一致。所要0.475秒。
- exceptional_geometry は (3) の全項と最後の六括弧を独立に手展開し受理した。
  generic_structure も最後の shear と実族・閉性の論理を独立に確認した。
  これらは別実装による再実行ではない。

入力 SHA256 と実行時の `producer_bytes` は証明書内に保存した。
旧二次曲線・三次曖昧さ商の診断証明書は不変であり、その生成時コードは各証明書の
`producer_bytes` から復元できる。旧 job を現在の拡張コードへ再 bind しない。
それらの未決定診断は今回の全支持証明の依存ではない。

fresh checkout での検証入口は、既存の出力を保護した一時コピーで
上記 job の gate を通した後 `--aff32-limit` を実行すること。
既存出力がある場所では runner は上書きを拒否する。
