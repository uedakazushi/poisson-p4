# generic B18のactual odd一次項とNの部分射影像：独立検証

`AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63` を受理する。対象は固定されたgeneric B18代表点でのactual CA/CB/N sourceと、共通29 target座標への **Nの形式的な部分射影像** である。Nの全target像やactual交叉Dの完成環を同定したものではない。

## 物理sourceと独立な一次項計算

別実装のraw CA/CB式、divergenceによるcanonical化、全150係数抽出、物理GL Lie作用から二つの全finite core mapを再構成した。CBの μ f の混合項も含む。十九のGL方向、九つのfiber方向、kernel旗方向のaugmented微分を再構成し、両sourceの29次元étale chartを確認した。保存されたproducerを新しい出力先へ一回実行し、全証明書がbyte一致することも確認した。

一次odd係数には、producerの二段階の射影と異なる一つのblock公式を使った。sourceの全150出力をFとし、中心を引いて必要な項を

\[
F=Js+Wp^2+P_s s p+Zp^3+\cdots
\]

とする。ここで (s) はkernel以外の28変数である。kernel方向の一次項が零なので、外側のGL方向とのmixed (p s) 項は零である。共通actual座標 (C=(Y_{28},r_{61})) を選び、

\[
M=[J\mid W],\qquad T=M_C^{-1},\qquad Q=1-MT\operatorname{pr}_C
\]

とおく。するとactual raw座標順のodd係数の一次行列は厳密に

\[
Q[P_s\mid Z]T
\]

であり、discriminantの微分はTの最終行である。この式は、(s,p^2) の一次変換を同時に行い、そのeven寄与を除く。展開で (s) に重み2、pに重み1を与えると、odd係数の一次部分に寄与する項は上記のmixed項とcubic項で尽くされる。全28変数に対するmixed係数と全150 cubic係数を直接抽出した。rank2基底のtraceを引く変更は、中心でodd係数が零なので、この一次部分を変えない。

この独立block公式の全150×29 entryは、producerがimplicit二次補正とcubic補正から得た行列に一致した。CA、CBともactual共通29方向のrankは7であり、二つのdiscriminant微分は独立である。選択29行自身のodd係数は零になることも確認した。これはfull odd idealの生成数を7に制限する結論ではない。

**座標に関する明確化**：保存されたtangent行列のrank6は、adjusted座標 \(\eta=r_{61}-LY\) を固定したときのY28方向である。actual raw \(r_{61}\) を固定した最初のY28列のrankは、CA、CBとも7である。実際、後者の列に \(b_{r_{61}}L\) を加えると前者のrank6行列になる。author noteの「Y28方向rank6」はこのadjusted \(\eta\) 固定の意味に限定して受理する。二つを同じraw部分空間とは扱わない。この明確化はNの以下の部分射影像の証明を変更しない。著者の追記 `AMBIENT_B18_ODD_FIRSTJET_COORDINATE_SCOPE_20261103_v63.md` をhash固定し、CA/CB/Nのactual行列をintrinsicな \(\ker d\nu\) へ制限したrankが全て6であることも別に確認した。Nの全28方向のodd rankは7である。

## Nの同じactual Y28を用いた底

CA sourceで (k_1=0) を課し、moving kernel旗を全次数で逆作用した後の全150係数が、固定旗N19空間に入ることを確認した。binaryの対応速度も零であり、この制約はjetだけの条件ではない。Nのaugmented source微分のrankは28で、proper N sourceの全完成chartを与える。

Nに対してはkernel以外の27方向とそのp²係数から28×28のblockを作り、そのactual target行を共通Y28そのものの順序に取った。これにより、先のN独自のregular底と同じactual座標集合であることを、単なる抽象同型ではなく確認した。独立block公式から得たextra row61のodd係数bとdiscriminant νについて、

\[
db=-dY_{16}-6dY_{65}+4dY_{105}+dY_{112}-4dY_{117},
\qquad \partial_{Y_{27}}\nu=-1/3,\quad
\partial_{Y_{27}}b=0
\]

を得た。従って (d\nu,db) は独立であり、scalar chartの (dY_{106}) と合わせたrankは3である。中心で (Y_{106}=1) なので、その射影chartを取っても二つの微分の独立性は残る。

## 全次数の部分射影像を得る代数的理由

既存の独立受理済み有限性と長さ2のparameter fiberから、Nの全完成normalizerは、actual底

\[
R=\mathbf C[[Y_{28}]],\qquad B=R[t]/(t^2-\nu)
\]

上のrank2自由algebraである。中心での座標は全て中心化している。extra actual targetの像は、一意に (a+bt) と書ける。この全係数はfinite sourceのimplicit消去とrank2基底による除法で定まり、一次項から全係数を補間しているのではない。

(z=r_{61}-a) と置くと、写像 (R[[z]]\to B)、(z\mapsto bt) のkernelは厳密に

\[
(z^2-\nu b^2)
\]

である。実際、monic二次式によるWeierstrass除法後の余りは (r_0+r_1z) であり、その像が零なら、BのR基底 (1,t) により (r_0=0,r_1b=0) となる。Rは整域、bは上記の非零一次項を持つので (r_1=0) である。従って追加kernelはない。

微分の独立性から ν,bをRの二つのregular parameterに選べる。全formal部分射影像はconeで

\[
\mathbf C[[g_1,\ldots,g_{26},\nu,b,z]]/(z^2-\nu b^2),
\]

projective chartではgが25個となる。νはparameterで非平方なのでこれは整域であり、bが非零なのでそのfraction fieldはBと同じである。Bは正則かつ有限なため、この部分射影像のnormalizationでもある。従ってordinary pinch pointと滑らか因子の積という結論は、全次数のrank2 algebraから従い、有限jetのみからの型の推測ではない。

全N targetにはさらにodd情報を持つ座標がある。CA/CBの相対image kernelの和、Dのnilradical、primary型、depthをこの部分射影から推論しない。形式的な完成射影像の結論を、全global射影閉包の結論へ拡張しない。

検証codeは `code/ambient_B18_odd_projection_review_20261103_v63.py`、receiptは `reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json` である。producerの再実行先には未使用のwork directoryを指定する。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/ambient_B18_odd_projection_review_20261103_v63.py \
  --replay work/B18-odd-projection-producer-replay \
  --output work/B18-odd-projection-independent-review.json
```
