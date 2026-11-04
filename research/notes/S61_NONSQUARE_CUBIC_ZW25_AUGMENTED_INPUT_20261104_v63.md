# projectedZW 25×31 とfullHpsi列：実25×32入力

2026-11-04。これはfull cubic56×90とは異なる射影入力である。
exact source chainは次の通り。

    full source → Z=W=0 → quad20+genuine16 の36×40
      → L座標、q1/q3をL²,Lx0へ変換
      → unit monomial9 source/target pivotとquad1零行を除く
      → quad10+genuine16 の26×31
      → literal duplicate G6=G1を除く →25×31.

31 targetはrestricted pure35から最初の九 monomial列を除いた26列とH*(L,x0,y1,y2,y3)の五列。
heavy5のSchurはまだ行っていない。
G6とG1の等式はZW restrictionと九 monomial projection後のものであり、full sourceではない。
実full40の差はL²*linearの五係数だけで、sourceの等しいclearing delta=−64000C²を使う。
full HZ/HW係数はZW restrictionで失われるので、full Cdeltaが非零であることと矛盾しない。

保存済みnonzero34minorは36×40rank34を示す。unit9消去でrank25となり、零行の削除ではrankを変えない。
従って25×31とそのaugmented25×32はgenericrank25。新たな点を評価したり、Delta10 determinantを展開したりしない。

## 追加psi列のliteral source

受理済みunitpairでCdeltaのfullH係数を(z_delta,w_delta)とする。
psi=−w_delta*(HZ coefficient)+z_delta*(HW coefficient)はCdeltaを零に送り、quad productsも零に送る。
元genuine行iのfullH係数二つを凍結済みraw57から直接読み、

    psi_i=lambda_i*(-w_delta*rawHZ_i+z_delta*rawHW_i),
    lambda_i=(-40)^3*C^(3−saved_C_cancel_i)

を構成する。L座標変更は五lightだけなのでZ/Wの係数は混ざらない。
G31にも使われた同じsource clearingを必ず掛ける。
quad10のpsiは零、G6とG1のpsiもliteralに一致する。
従って同じG6一行だけを除いたactual25×32入力が得られる。

rankaugmented=25なら、既知product41とglobalCdelta unitimageによりfullgenuine16は独立となり、
**fullcubicrank57の十分条件**になる。projectedgenuine rank14ならCdeltaを除いた残りkernel一行に
psiが非零であることでfullrankを救える可能性がある。
psiが零でもpureT27の像が非零ならfullrankはありうるので、augmentedrank条件は必要条件ではない。
projectedrank≤13ならpsi一列では全missingdirectionsを扱えず、pure directionsが別に必要。
選択したDelta34が零であることと、projectedrankが落ちることも同一視しない。

初段は既存25×31係数をコピーし、新たには16個のpsi scalar coefficientsだけを展開する。
input100k core terms・32MiBのguardを通るときだけ、別の承認済みONE known-unit Gaussianを実施できる。
このinput packet自体はpivot/rank/GBを実行しない。

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_cubic_zw25_augmented_input_20261104_v63.py --output work/marked-square-v63/cubic-zw25-input-REPLAY
```
