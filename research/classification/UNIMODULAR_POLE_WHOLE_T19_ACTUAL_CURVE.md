# h5残一軌道はT19：自然graded moment実曲線と固定符号正逆

2026-11-08、PRIMARY E_multiple_hyperplane。producer全実曲線証明、独立別作者査読前。
新ambient／point／minor／local jet無し。旧D101/M47非閉鎖を再開せず、真T19の固定whole source実曲線から帰属を証明する。

## 1. 一度だけ実行したuniversal曲線と限定的FALSEの正確な意味

`CLASSIFICATION_UNIMODULAR_POLE_T19_GRADED_MOMENT_CURVE_20261108.json`
は固定frame／全R,V／全inverse／det／true completed source／rawEuler／全150canonical Laurent／
strict10／元350／全126actualPfを一括に検証した。全exact checkはPASS、3.241401946秒／55,028KiB。
全150canonical係数の負e冪は空で、実曲線は全e=0まで正則である。
raw liftには純Eulerのe^-1項が残るため、raw lift自体の正則性は主張しない。
保存したcanonical sectionと全rawEuler差でactual classの正則性を証明する。

callerの最後の自動criterionは中心を直接aR+bBに入れることだけを許した。
保存した完全線形系はV=0を強制し、その限定criterionはFALSEとなった。
この報告／caller／jobは一byteも変更しない。
下記の**固定符号actual正逆**により中心は同じwhole h5の一実軌道に戻る。
新curve係数又はfree frame parameterを増やしたのではなく、既normalizerの許す定数GL4を手で適用する。

## 2. 全保存中心と二つのrank二tensor

h行零の六成分を順(xy,xz,xw,yz,yw,zw)で書く。

```
R=(hx,−x²,xz−h²,xy,−yz+hw,−xw+hy),
Rtilde=(hx,−x²,xz+h²,xy,−yz+hw,−xw−hy),
B=(0,0,0,−h²,−hy,xy+hz).
```

R=Jac(xy+hz,xw+hy+z²/2)、Rtilde=Jac(xy+hz,xw−hy+z²/2)（h一chart）、
B=Z1∧Z2である。
保存全150係数の中心は正確に

```
P0(Rpar,V)=−(V/2) Rtilde−Rpar B.                    (1)
```

全表は

|pair|保存中心|
|---|---|
|xy|−V hx/2|
|xz|V x²/2|
|xw|−V(xz+h²)/2|
|yz|Rpar h²−Vxy/2|
|yw|Rpar hy−Vhw/2+Vyz/2|
|zw|Vhy/2−Rpar hz−Rpar xy+Vxw/2|

であり、h行は全零。
一つのPfaffianや選択sliceの一致ではなく、保存全150basisの非零十三係数をそのまま再掲したものである。

## 3. 固定QQ符号frame

```
S:(h,x,y,z,w) -> (h,x,−y,−z,w), detS=1, S^-1=S.
```

全六式の直接置換で

```
S_*Rtilde=−R, S_*B=B.                                (2)
```

例えばRtildeのxwはold xz+h²が−xz+h²へ、zwは(−1)(−xw+h y)=xw−hyへ移り、
どちらも−Rの対応成分である。Bのyzは符号二個で不変、yw／zwも変数符号とtensor符号が相殺する。
従ってS_*P0=(V/2)R−Rpar B。
特に**実行前に予測した同じ定数(Rpar,V)=(-1,2)**を用いると、中心は厳密に

```
S_*P0=R+B=h²C+L+Q=P*.
```

追加の数学的点を選んだのではなく、既保存universal Rpar,V恒等式の一つのunit値である。
全Rpar Vunitについてもa=V/2,b=−Rpar unitであり、事前のactualdiag/tensor normalizerで同じP*へ戻る。
本稿の正結論には上記QQ定数だけで足りる。

## 4. 真T19からP*への全actual曲線

true completed T19のsource順(a,b,c,d,H)、period(r,v)、edges(t0,t1)を使う。
必要なcoupled correction Pac=2t0t1 bH/(r−v)を保持する。

```
r=−e³, v=2e²−e³, t0=−e², t1=−2e².                 (3)
```

period openは−r(r−v)=−2e5非零（eunit）。
全actual source→最終target frameはSと保存momentframeを合成したもの：

```
h=−2e³(b+c), x=e⁴(a−2b−4c), y=−e(2b+c+d),
z=e²(a+2b+2c), w=−a+2b+c+d+H.                       (4)
```

det=4e10、逆は

```
a=h/e³+z/e²,
b=x/(2e⁴)−3h/(2e³)−z/(2e²),
c=h/e³−x/(2e⁴)+z/(2e²),
d=−y/e+2h/e³−x/(2e⁴)+z/(2e²),
H=w+y/e+z/e²+h/e³.                                  (5)
```

全eunitでactualPGL、追加tensor scale一、二つの合成はidentity。
保存した全universal150canonical係数は負冪無しで、(3)と固定符号を適用しても正則である。
e0の全150中心は§3のP*。generic eは真completed T19のactual像である。

fullactual divisorは、(5)を代入した

```
4e10 * [-r(r−v) H(cd−b²)(ab+dH)]                     (6)
```

であり、負冪無しは保存全126係数／共変則で確認済み。
中心は−Rpar V h5/2=h5、P*の元signedEulerPfaffianに一致する。
rawEuler場は保存した全divergenceをSでpushするだけで、canonical sectionとの差を全項保持する。

## 5. 閉性によるwhole被覆

eunitでは(3)はperiod open内の真T19全族、(4)はactual frameだから、全π_eは名前付き閉成分T19に入る。
全150canonical係数はe多項式、π0=P*非零。したがってP*∈T19。
T19はclosed PGL不変なので、既one-orbit normalizerの全unit逆と四translation／最高GL4を戻して、
元x²y²alignedzero全sourceの全点がT19に入る。

`UNIMODULAR_H5_WHOLE_ONE_ORBIT_SEPARATE_REVIEW.md` が証明した
全unimodularconstant-nonzero-Pf C6/L24/Q60の「既知閉和＋一軌道」へ合成すると、
**全A=c=0のrank四h5源が既知閉成分の和に入る。**
最高Q零／factor／nil／non-split、theta非primitive、全低次unit／zero補集合を除外し直さない。
Gは既に閉じており再活性化しない。新component増分零。

本稿はproducer証明であり、別作者literal全150/rawEuler/126同曲線の独立再生は別receiptへ分ける。
全classificationschemeのsmoothness、P*の全局所枝／nonreduced厚みを主張しない。
