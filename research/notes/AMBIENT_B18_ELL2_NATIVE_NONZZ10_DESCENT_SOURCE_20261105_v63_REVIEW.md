# native nonzz10 の実 source 降下：独立検証

独立 SOURCE 検証は PASS した。報告 SHA256 は
`79a6ba0ff2ccc0190e8ac4d8a46515514349a75fa62428e28713a12a91a99ff6`。
元の64個の IA×IB 二次初期式、33基底上の10 functional、実660単項式の
全6,600個の NF 恒等式を確認した。kappa、G8、honest quotient、応答値、
新しい jet、upper 値は読んでいない。

共通の実座標は zA0,…,zA6,rA,ell、ell=5(sA−rA) である。
元 IA の八生成元は CA の七座標と raw f63=node7790、IB は CB の七座標と
raw f63=node28539である。両 raw f63 がそれぞれ自身の中心付き F9 chart で
定数・一次を持たないことを、既存 source2 の literal coefficient から確認した。
source2_monomials は指数ベクトルでなく variable-slot multiset なので、
その記法を独立に九指数へ変換した。中心付きの実共通座標への pullback はこの次数を保存する。

実 CB 一次変換の最初の七行に rA,sA 成分はない。
従って IA の一次式は zi または零、IB の一次式は A0 zi または零となる。
全64積の二次初期式を元の source ID の順に再構成すると、全て zz のみである。
これは正規化後の模型だけの議論ではない。任意の形式係数による積の二次項に
寄与するのはその定数だけなので、十個の nonzz 二次係数 functional は
実 product ideal 全体を消す。ell²と m³も消す。普通の Taylor 係数を取り、
二階微分の factorial は導入しない。

十 functional の対象は z0r,z0ell,z1r,z1ell,z4r,z4ell,z5r,z5ell,r²,r ell であり、
実 raw33 基底上の selector は順に
22,23,20,21,14,15,12,13,24,25 である。
受理済みの source equality と m⁵=0 の証明により、この raw33 は
実 R/(IAIB,ell²) の商である。

非斉次 carry を独立に点検するため、保存済み660 NF と functional の積だけではなく、
元 raw33 の九座標の乗法作用から660 NFを全て再計算した。
九作用の共通分母 D=11828436480 を取り、整数行列で一単項式ずつ作用させる。
出力の次数を切り捨てず、非斉次項をそのまま保持する。次数四以下かつ ell 指数二未満の
全660入力について、再計算結果が元 NF と一致し、全6,600個の functional 値が
通常係数の Kronecker delta に一致した。

事前 guard は積・NF検査より前に保存した。整数作用の最大列絶対値和 L から、
四段の分子は L⁴（150 bit）、分母は D⁴（134 bit）以下となる。
各積・部分和も同じ正値上界で管理し、保存済み NF との比較も別の事前上界で管理した。
実際の比較は134 bit以内だった。

原報告の prior1 は selector–NF 検査の範囲である。64積の係数は CA 側が単位座標のため
CB 一次係数のコピーで、その別の上界は7 bitである。原事前 phase の記述と
不変の HEIGHT clarification を照合し、この二つの範囲を混同していない。
元の slot-multiset reader の停止履歴も保存された hash と一致する。
原資料のバイト列や予算は変更していない。

OLD Q9 の C0 基底末尾は gA=node505938であり、raw f63=node7790とは区別した。
同じ実 F9 順序・ell chart・既存 current M0 transport を照合し、
M1 や比較 gauge を零とは置いていない。
将来 native canonical kappa/G8 に functional を成分ごとに適用するには、
同じ branch・monic・free4・overlap・資源の別 SOURCE gate が必要である。
raw33 は非単位である g の除算を許可せず、honest g/h 商の像も自動には与えない。

再現 checker:

```sh
.venv/bin/python code/ambient_B18_ell2_native_nonzz10_descent_source_review_20261105_v63.py
```

記録済み実行は0.108秒。source所有数上界66,057、将来状態・ラベル・予備領域込み
118,557で、150,000 slots、60秒、2 GiB、2048 bit、RAW32 MiB以内である。
凍結された同じ出力先への再実行を拒否する。別作業コピーで再検証すること。
