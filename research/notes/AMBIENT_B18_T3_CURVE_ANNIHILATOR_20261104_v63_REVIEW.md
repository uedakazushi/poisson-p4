# T3のactual曲線annihilatorと一parameter判定：独立監査

2026-11-04。**I0_CA*T3=0、multiplicity3のCM曲線、明示parameter ell、閉fiber長さ≤36、
一colonによるT3有限台torsionの消失判定を受理する。** colon自体やT3の零非零、全Dのdepthは未決。
authorやsourceの再実行、新jet、点、GBは無い。

## 96恒等式は元のactual quotient内で成立する

full186入力の独立reviewと保存matrixを使う。Q[3]のlower部分は180生成元×H十一列で1,980列、
その後は六cubic×C0九列。従って0-basedで
`column=1980+9alpha+(i+1)` は実際にc_alpha*(f_i*tau_A)を表す。
top列2alpha+eはc_alpha*(j_e*tau_A)。

独立な非再帰的product normalizerで、全96個×二行について
`f_i*top_(alpha,e)=j_e*Q3_column` を可換積のfactor multisetとして確認する。
非乗法nodeはopaqueなactual函数として扱い、有限jetや未知係数を零にしていない。
全pointer、f_i、j_e、両行signatureを保存recordへbindする。
十二top列がT3を生成する既受理の結果から、I0_CAがT3全体を消す。

これはambient R0²をI0で割ってmatrix imageを取る操作ではない。
例えばJ=(x),Q=(x²)ならJ/Q≅kは非零だが、ambient columnsはmod xで双方零になる。
Torによりこのような像の取り直しは誤るので、quotient-moduleの恒等式のまま使う。

## Actual CA曲線とmultiplicity

元source DAGの独立Fraction ADからCA Hensel48のimplicit derivativeを求め、
実10×9 odd Jacobian全成分を再計算する。
detH48=-84934656/625、選択(f18,f19,f28,f29,f60,f67,f72,j46,j95)のdet=128/3。
f63の一次部分は零。この九本が同じphysical R0の正則座標(z7,r,s)となる。
一つのpointed57 inverseの全式・bindingもread-onlyに再構成し、
g_A=f63|z0が元のactual sourceへ結び付いていることを確認する。

coordinate telescopingによりI0=(z7,g_A)。g_A=0ならI_CA=I0+wJは九本で生成され、
受理済みactual μ(I_CA)=10に反する。従ってC_A=k[[r,s]]/(g_A)はCM一次元曲線。

actual conductor quotientを同じG19で切った環は
`k[[w,r,s]]/(g_A,wr,ws)=k[[w]]×_k C_A`。
最大idealの正冪は二因子の直和なのでHilbert–Samuel multiplicityは加法的。
共通群productがconductorにも作用するため、受理済みmultiplicity4はG19切断後も保たれる。
従って4=1+e(C_A)、e(C_A)=ord(g_A)=3。
旧analytic模型のwやsr²を新しいactual座標と同一視していない。

## Parameterのphysical tangent binding

checkerはauthorのSymPy処理とは別にFraction Gaussianで、旧29×22 conductor tangent frameから
同じphysical G19切断のSchur項を引く。
全29行で群19列が消えること、旧UVのw速度が零、七zの微分が零、旧Wのw速度が非零を確認する。
そのUVからactual(r,s)へのblockは

```
[[4/25,-4/25],[9/25,1/25]], det=8/125。
```

従ってexact函数ell=5(s-r)の一次部分は旧U+V。
受理済み**全associated-graded ideal** `(WU,WV,(12/125)U(U-V/3)^2)` を使う。
新しいwはUV上で一次零、(r,s)のUV blockは可逆なので、新曲線の三次初期形は
この旧binary cubicへ非零scaleで対応する。高次座標変更の加速度を同一視する必要はない。
full initial idealのdegree3でWを落とした主idealからこの主張が従うので、古い一点の
cubic係数だけを新しい曲線へ外挿していない。

ell=0の旧方向(U,V)=(1,-1)における値は64/375≠0。
よってellは曲線tangentを含まず、ord(g_A)=3からlength C_A/(ell)=3。
C_AはCM一次元なのでellはnonzerodivisorで、Weierstrassによりk[[ell]]上free rank3でもある。
exact ellは同じfull186 source namespaceの二nodeの線型結合として保存されている。

## 一colonの消失判定と長さ36の限界

I0*T3=0かつrad(I0,ell)=m0なので
`H0_m0(T3)=H0_(ell)(T3)`。
ell-power torsionが非零なら最小の消失冪の一つ前を掛けて非零ell-kernelを得る。
従ってその消失は
`((Q[3]:ell)∩Jinf)/Q[3]=0` と同値。
ell-kernelの次元をm0-socleの次元と同一視してはいない。
colonは元のactual R0² submodulesで行う。

十二generatorと曲線ell-fiber長さ3から、T3/ell*T3の長さは高々36。
このfiber上界はT3自体の長さ、torsion exponent、jet次数の上界ではない。
一生成のk[[ell]]/(ell^N)でもfiber長さ1に対しNは任意に大きい。

今回colonは未計算。T3の有限台torsionが見つかればcone depth19となるが、
このcolonが零でもT0,T1,T2を省略してdepth20と断定できない。

再現:

```
python code/ambient_B18_T3_curve_annihilator_review_20261104_v63.py \
  --output /tmp/B18-T3-curve-independent-review
```
