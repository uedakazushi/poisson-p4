# Actual CBの四g商の原点値：独立microjet監査

2026-11-04。**四つの既知正則商の原点値が厳密に0であることを受理する。**
商全体の零、必要な高次jet、T3所属・colon・depthの結論は含まれない。
凍結済みfull DAGのNone値はそのまま保つ。

## 同じ57系を独立に再構築する

actual CB sourceの48 Hensel式と、actual F9=(zero7,j46,j95)の九式を使う。
未知順は元のsource48係数とphysical X9の計57。
targetは(0,...,0,r,s)で、別のsource点やungauged系へ移していない。
独立な全DAG Fraction ADから57×57 Jacobianを構成し、保存matrix全成分と一致させる。
det=-3456/625を独立Gaussianで確認し、保存J逆行列の左右両積も単位行列へ照合する。

元41315 nodeのうち、全57式、f63、四bbarに必要な9978 ancestorを独立に求め、
保存した全node IDsと一致させる。
microjetはauthorの固定10配列/recursive inverseと別に、非零monomialのみの辞書で表す。
単元逆元は定数項で割った後の有限幾何級数1-h+h²-h³として構成する。
QQ[r,s]/(r,s)^4内で次数1、2、3を固定J逆行列で順に補正する。
全homogeneous補正、全57未知数の10係数、元57式の全10係数を独立に一致させる。
authorの実行や高次solverの再実行ではない。

## cubicと正則商の定数項

独立に得たf63|z0の三次初項は

```
-24r³-48r²s+18rs²+54s³ = -6(r-s)(2r+3s)²。
```

lower degreeは零で、r³係数-24が非零。
四zero physical rows64/68/73/74のbbarは次数0〜3の全係数が零。
それらのfull DAG g商nodeへの対応もlabel/原source rootで照合する。

重要な入力は、前のFULL_RELATIVE_INPUTレビューが全次数で
`bbar=g*h`、hが正則algebraic power seriesと証明済みなことである。
したがって最低次成分は `bbar_3=g_3*h(0)`。
g3≠0とbbar3=0から四つともh(0)=0が厳密に決まる。
有限jetから正則性や全次数membershipを推測していない。

今回の値は独立なoverride receiptへ保存するだけで、元のfull DAGの
`not evaluated; no unit use`という歴史的metadataは変更しない。
gがこの三次多項式に等しいこと、hが恒等零であることは言わない。

## 次のactual-colon計算への条件

四定数が決まったことで、full186有限回路の原点値は派生receiptで正当に伝播できる。
そのresidue/unit-pivot preflightと、同じintrinsic T3に対するCB側annihilatorの検査は
manuscript担当と調整し、重複実行していない。
literal二行行列のresidue rankだけではQ[3]への所属やT3のsocleは分からない。
I0によるambient reductionにもTorがあるため、actual quotient/presentationを保つ必要がある。
次の消失判定は既受理の `((Q[3]:ell)∩Jinf)/Q[3]` であり、その計算は未完。

再現:

```
python code/ambient_B18_CB_g_division_origins_review_20261104_v63.py \
  --output /tmp/B18-CB-g-origins-independent-review
```
