# nil31短translation・H+既約quarticの全MPS被覆

2026-11-07。全源の有限化と実際のPoisson族による被覆証明。
主課題Eの型4 `Z=y∂x+z∂y+∂w`。新成分を追加しない。

## 全coreと全quartic補集合

既済の全九core枝を用いる。記号はc=n17,a=n11,b=n16,t=n15,j=n10,k=n9。
c=t=k=j=0の連続二枝では、最初のlower式が
`(3a+1)n8/6=0`（b=1/3）または `(1−3b)n8/9=0`（a=−b）。
n8=0のもとで最初のmagnetic式は `(a+1/4)n7=0`、最高Pfは `a n7 z4`。
従ってquartic条件は三組 `(a,b)=(-1/3,1/3),(-1/4,1/3),(-1/4,1/4)` のみに残る。
これと非零c五枝・非零t一枝・非零k一枝で全10候補である。

`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_SOURCES_20261107.json` は、全10候補でlower/magneticの
二段線形系を一様なmonomial pivotで解き、全残compatibilityを保持した。
非零cの第四枝以外はそれぞれ z、z²、3cz+1、15cz−1、z+3t、highestPf0、z²、z、z4
により、全parameterで既約quarticを持たない。
唯一生き残る完全源は

    c!=0, a=-5/12, b=-1/12, t=-1/(72c),
    j=7/(72c), k=-1/(216c²),
    free=(c,n12,n13,n6,n7,n8).

全Jacobiが恒等零で、parameter空間はGm×A5。
Pfはx一次で係数 `c(z−1/(6c))³`、その唯一の零点z=1/(6c)での残りは

    −(12n12c+2n13+6cy−w)²/(3456c³).

w²係数が非零なのでprimitiveなx一次式で、全parameterで絶対既約。
この全源還元は `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_SUPPORT_20261107.json` に保存済み。

## 実際の有限Poisson族

c=1、他の五freeparameter=0の点の整数432倍をPとする。
h-Casimir代表の六括弧は

    Pxy=−w²+3wy−108xz+36hx−36y²,
    Pxz=3wz−2hw−180yz+39hy,
    Pxw=−3w²−180wy−1296xz+324hx+864y²,
    Pyz=−2h²+42hz−180z²,
    Pyw=−180wz−6hw+432yz+144hy,
    Pzw=432z²−36hz−6h².

実群/gauge/scaleの47方向に横断するRと補正Tを

    Rxy=0,
    Rxz=126xy+6xw,
    Rxw=648x²,
    Ryz=−132hx+792xz−18y²+6yw+w²,
    Ryw=−2376xy+504xw,
    Rzw=684hx−4104xz+864y²−126yw−12w²,

    Tyz=−504x², Tzw=3024x², その他Tij=0

とする。全350式で **P(e)=P+eR+e²T が全eについてPoisson** と検証した。
これはformal jetではなく、表示した多項式族そのものの恒等式である。

発見に用いたC*重みは(h,x,y,z,w)=(0,2,1,0,1)。P,R,Tのtensor characterは−1,1,3。
`diag(1,s²,s,1,s)_*P(e)=s P(e/s²)` なので、非零eは全て実PGL/scalar同値となる。
最終成分認識には高次のformalモデル完全性を仮定する必要はない。

## 全周囲唯一性とMPSへの帰属

全Dは350行150列でrank102、vertical/実群/scalarはrank47、さらにRを加えるとrank48。
したがって実Poisson族の射影次元下界は23、全接次元上界も125−102=23。
Pは滑らかな23次元点で、通る周囲成分は唯一となる。

P(e)はh-Casimirを保ち、全homogeneous modular cofactorは

    Zx=432y,
    Zy=432z−1296e x,
    Zz=e(216y+36w),
    Zw=432h−5184e x,

であり、特性多項式は

    u(u²+93312e)(u²+373248e).

従って一般eはregular semisimple、affineスペクトル比は(-2,-1,1,2)。
これは既済nonzero-semisimple入口の型12で、全H+既約quartic支持が有限閉MPSunionに入る。
e=0のquarticが既約なので、quartic既約性とhを因子に持たない条件は一般eで保たれる。
従って一般P(e)は既知MPSunionに入り、閉性によりP(0)も入る。
P(0)を通る成分が唯一であることから、その23次元成分は既知MPS成分である。
最後に完全六源Gm×A5は既約でこの点を含むので、全像とPGL飽和は同じ閉成分に入る。
これで型4の全H+既約quartic支持が既知MPSへ覆われる。

## 証明書と検証状態

全10候補と全線形解は
`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_SOURCES_20261107.json`、
九候補の全排除と唯一六源の全既約性は
`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_SUPPORT_20261107.json` に保存した。
後者はc以外の反転を持たないGm×A5の全源であり、連続coreの補集合も含めて全支持を覆う。

最終の実Poisson族・全rank・cofactor・既知MPS認識は
`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_RECOGNITION_20261107.json`。
reusable runner `code/classification_hyperplane_quartic.py` の
`--nil31-short-quartic-recognition`、有限契約
`jobs/hyperplane_quartic_nil31_short_quartic_recognition.json` で再現する。
最終報告はwork内の探索結果を参照せず、全350式の実curve恒等式と非零102/48小行列式を直接検査する。
既存semisimple型12の全支持MPS被覆は、その数学ノートと確定証明書のhashへ束縛した外部依存である。

`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_REVIEW_V2_20261107.json` は
別のFraction/Laurent実装で、元の全90係数・110線形行、全10候補の二段unit pivot、
残compatibility n13*n8、九つの普遍Pf排除、唯一六源の全係数とGauss既約性を再構成した。
さらに整数432P、表示実curveの全350式と全10raw Jacobi、D102・vertical24・M47・(M,R)48、
全cofactor/特性多項式/C*恒等式、cofactor軌道13を独立PASSした。
全源の射影飽和は6+13+4−1≤22、一方P0の実orbitは22なので、既約な全源飽和閉包はP0のorbitclosureそのものである。
この事実とMPS有限unionの閉PGL安定性だけでも全源被覆が従う。
滑らかな23次元唯一成分という上記の補強も全rankから独立確認済みである。
旧二次/三次support診断は結論を出さず、最終命題には表示した実多項式Poisson族を使用する。
高次jetやnilpotence指数の決定を分類の前提にしていない。
