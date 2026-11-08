# 二simple H+cubicの残m3二源：T37の一つの全rational plane入口

2026-11-08、PRIMARY E、`/root/quartic_sources`。純手監査。新演算、kernel、operator、point、profile、GB、jetは無い。対象はnode-flex/tangent-through-flexとcusp-flex/secant-through-flexの**全六磁気源**であり、既正則shear後の全lambda表示は `TWO_LINES_PLANE_CUBIC_THREE_MAGNETIC_QUOTIENTS.md` を使う。第三のm4源は別担当である。

## 1. 使える既存literal入力と歴史的未確認情報

固定した真completed T37 の全period/smoothing一のtensorについて、`CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_OPEN_LOCUS_20261107.json` はcell012の実plane行列、全六Lagrangian式、全restricted actual quinticを保持する。

```
oldX0=A, oldX1=a0*A+a1*B+a2*C,
oldX2=a3*A+a4*B+a5*C, oldX3=B, oldX4=C.
```

保存された最後のGB式には `(a5+6*a4²)*B(a4,a5)` という二factorがある。ここでB側の全actual sourceがLQ²であることは `T37_LAGRANGIAN_CURVE_CLOSURE.md` で既に別に閉じており、再開しない。

過去chatにあった「rational枝のsmooth conic+cubic、conic determinant324e³」という記述には、担当者のread-only確認で独立immutable報告が見つからなかった。その数値・因子型を証明の前提にしない。本稿が使うのは上記OPEN_LOCUSの実literal入力だけである。

## 2. Gm一曲線の全actual planeと逆

u≠0とし、次の六係数を取る。

```
(a0,a1,a2,a3,a4,a5)
  =(u/2, -1/(216u²), 1/(36u), -1/(72u), u, -6u²).       (1)
```

全original六式へ直接代入すると、順に次の零恒等式を得る。各式の保存時の非零constant倍を除いているだけであり、必要条件の一部を選んだものではない。

```
48*a0*(a2+a3)-24*a1*a5+24*a2*u-48*a3*u-1
 =1/3-2/3+2/3+2/3-1=0,
6*a0*u-3*a1*a2-6*a1*a3+3*u²+a5
 =3u²+1/(2592u³)-1/(2592u³)-3u²=0,
6*a0*a5+a1+3*a2²-6*a2*a3-3*u*a5
 =-18u³-1/(216u²)+1/(432u²)+1/(432u²)+18u³=0,
3*a1*u-a2-3*a3=-1/(72u)-1/(36u)+1/(24u)=0,
18*a1*a5-54*a2*u+1=1/2-3/2+1=0,
3*a0+3*a2*a5-u=3u/2-u/2-u=0.
```

逆はu=a4であり、(1)はsource上の全Gmを保持する。v=a5=−6u²にu=0を足すと、元第五式は1=0になるので、このcell内にはu=0の点は無い。Grassmannianの別cellでのprojective境界は別問題として保持する。

この一つの閉曲線の支持が因子 `a5+6*a4²=0` の全持上げsourceを尽くすとは本稿だけでは主張しない。Bとの交差に追加liftがあるか、他cellのboundaryがあるかを捨てない。今回必要なら調べるのは、あくまで(1)という明示的な**一つの全parameter源**である。

## 3. 最弱の有限必要義務とその限界

このsourceを既literal restricted quinticへ戻し、全21係数、factor積、degree二factorが得られる場合のdeterminantと全退化parameterを一度確認すれば、二simple線と既約cubicのmarked配置をこの一曲線が供給するか判定できる。有限contract・明示GO前にはこの計算を開始しない。新GBや新Lagrangian kernelは不要である。

generic factorがsmooth conic+cubicであっても、それだけで対象の二simple H+cubic源へ所属を主張できない。conicが分裂する全parameterとそのprojective境界を調べ、対象のnode-flex又はcusp-flex marked配置を実際に得る必要がある。さらに、zero-magnetic seedの一致だけでは不十分で、元四次元translation像の二次元cokernelを全て保持するactual moving-period/plane曲線が必要になる。

既T36 mainのuniform normal-cone入口は、binary Q2とcubic C3が同じsingular pointを持つproper閉条件を満たす。対象の第一二配置はこのmarked条件を満たさないので、その全branchの単純な再利用は閉じない。これはT36全成分の陰性ではなく、非一様moving actual frameや別period曲線を排除していない。

新 `TWO_LINES_PLANE_CUBIC_ZERO_LOCUS_T38_EXCLUSION.md` のgeneric-lambda零scheme論証はproducer純証明として別保存している。special lambdaのzero-planeを除外せず、独立受理前に既知親全分類と混同しない。現時点で第一二whole源のpositive parentは未決である。
