# nil31 short r=2/3 全源：真T19 moving-H曲線の別作者査読

2026-11-08、PRIMARY E、別作者 `/root/quartic_sources`。純証明 a087bd7818fc8b3d0ff57edc6d5c54ffecfa4800d3b9de856cd210b08093901f を手で監査し、明示GO・fresh gate後、同一 universal 曲線だけを一度厳密に再構成した。新しい点、周囲微分、operator、jet、period profile は計算しない。共有コード、旧証明書、有限契約の作成時bytesは変更していない。

## 1. 全源への矢印と零補集合

`MODULAR_NIL31_SHORT_SUPPORT.md` §7 の実affine逆は元六自由変数を全て保持して二変数(t,d)へ送る。rank四の条件は t≠0 であり、d は任意である。今回の中心は t0=1/27,d0=−sigma/108 を全sigmaについて直接与える。

既証明の実対角変更とbivector scalarによる作用

`(t0,d0) -> (mu*nu*t0,mu²*nu*d0)`

で mu=1,nu=27t,sigma=−4d/t と選べば、任意の t≠0,d を得る。sigmaもdも反転しない。従ってd=0を一般点だけの推定で補う必要はない。元六変数の全逆とこのtorusで全源へ戻る。

## 2. actual frameと全六中心の手監査

真completed T19 の coupled entry は `2*t0*u0*bH/(r−v)` であり、t0=(r−v)/2,u0=r−v によって `(r−v)bH` になる。r=e²,v=3e²+e⁴sigma の period開は `e⁴(2+e²sigma)`。各sigmaに対してgeneric eで非空である。

source-from-intermediate frame(3)を直接掛けると q3=hZ−e²W²,q4=hY−ZW+e²XW、親hyperplaneは H=Y+X−Z/e² となる。両quadricを保存するGa(lambda=1/e,mu=0)の逆向きと一致し、detはe²。forward式(4)では特に Y=H−a/e であり、符号をforward flowと取り違えていない。

chart六式(6)を手で追った。PXYのq除算では、(Z−qX)Kのq一次係数とPYZのq一次係数を引き、最後の `3Z(Y−ZW)` を加えると、ZWZの項が相殺して

`PXY|e0=2XZ+4YZ−6WY+sigma*Z²/2`

となる。PXWの第二のq除算も正確で、残五式と合わせて(7)を得る。Euler補正L=(0,−4Z+4W,2Z,0,2h)を引くとZ行が零になり、全六(8)へ移る。固定GL(9)の全逆とdet1/162を確認し、全六pushは(10)、すなわち t=1/27,d=−sigma/108 と一致する。特にPxzのnz係数1/27、Pxwのsigma係数−1/108、Pxyのsigma*nz係数−1/108を個別に確認した。

一般eで旧hをCasimirとはしない。chart restrictionのinjective linear imageとGL-equivariant canonical divergence-zero sectionを用いる議論は有効である。次の厳密readbackではこの議論だけに依存せず、全150係数自体に負e冪が無いことも確認した。

## 3. 同一universal曲線の厳密再構成

有限契約 `jobs/nil31_short_two_thirds_t19_whole_confluence.json` の全入力SHAとfresh gateがPASSした後、標準ライブラリFractionと既独立LaurentRingだけを使い、QQ[e,e^-1,sigma]上で一度実行した。CASとproducerコードのimportは無い。制限60秒/1GiB、実測0.158381900秒、最大RSS13,568KiB、全PASS。

証明書 `reports/CLASSIFICATION_NIL31_SHORT_TWO_THIRDS_T19_WHOLE_CONFLUENCE_20261108.json` に、実行inline全bytes、真source全十entry、period/correction掛戻し、Ga両inverseと両quadric保存、全frame正逆、全150 canonical curve、rawEuler全forward/inverse、chart全六と二つのq厳密割算、中心gaugeと固定GL、全126 actual Pfaffian、全十strict/350projective Jacobiを保存した。

全canonical150は座標について斉次二次、eについて負冪が無く、最大e次数は六である。最大二次等への不当な打切りは行っていない。全intermediate actualDは

`(2+e²sigma)*(e²(X+Y)−Z)*(hZ−e²W²)*(hY−ZW+e²XW)`

そのものであり、det変換とrawEuler不変性を別々に掛け戻した。固定GL後の中心は

`−n²*z*(n*y−z*w)/27`。

全canonical中心は目標n-Casimir六式から `E wedge N/5` を引いたものと全150で一致する。ここで N=y∂x+z∂y+n∂w。従ってstrict liftと目標Casimir表示を混同していない。

## 4. 受理範囲

全sigmaについて、正則な非零中心とgeneric eの真T19所属が得られた。named閉T19の閉性、§1の全torus逆、元六変数の全normalizer逆を合成し、**nil31 short r=2/3 のrank四全源は点ごとにT19に属する**。

これは全源の陽性既知親認識であり、新成分を数えない。旧二次/三次局所診断、固定一様blowupの不成立、旧frontierのbytesは変更しない。それらから成分性を推定したのではなく、別のmoving actual Ga/period曲線により閉じた。他の未解決divisor型や全分類の完了は主張しない。
