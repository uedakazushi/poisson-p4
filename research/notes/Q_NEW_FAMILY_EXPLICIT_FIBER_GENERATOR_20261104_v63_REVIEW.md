# Actual fiberの明示生成元：deckと同じ関数体の独立監査

2026-11-04。**k(Y16)=k(B15)(w_fiber)、Jのactual fiber上degree2を受理する。** k=QQ,C。
一般関数体の主張であり、全点・全special fiber・全域finite morphismを分類しない。
旧存在証明のLüroth生成元を、今回具体化したものである。

## Source deckはparameter同定やβ kernelと別である

独立review済みの全tensor Aut Eとsigma sign+1を使う。
checkerはexpression cancelのauthorとは別のsparse rational field演算で、保存行列のX,W成分と
spatial determinantを照合する。z依存成分が無く、modular GaのX,W rowが零であることも確認する。
Eはrf,sf両符号を反転し、sigmaはrfだけ反転する。
sigmaの実scalarとPfaffian則は `v'=-32u-v,c'=-cv²/v'²,B'=Bv³/v'³` を与える。

F5の五定数の定義へ両変換を直接入れ、全て恒等的に保存されることを独立確認する。
これらはfull projective tensorの既証有限同値であり、βだけの同値ではない。
正規化したfull tensor間の差はSp4 frameだが、同じ定数から同じβ sectionを得る。
その幾何的stabilizerはidentityなので一意なactual frame hも固定され、full tensor自身が固定される。

このsource deck V4は、normal form parameter(u,v)上の別のV4と同一ではない。
full sourceの次数8を数える際にも余分なfactorを掛けない。

## 固定体の計算とactual fieldとの等号

A=256、D=Kt²-K+1、N=Kt-K+1として
`rf=A(D-2N)/D,sf=A(D-2N)/(D-2tN)`。
二つのMöbius変換 `N/[K(t-1)]`、`(K-1)/(Kt)` のinvolution性・可換性・符号作用を
独立なsparse rational演算で照合する。異なる非自明な二involutionなので群は忠実V4。
rf²は固定され、tについてdegree4（約分後のnumerator/denominatorで確認）だから
`F5(t)^V4=F5(rf²)`。

R=rf²,S=sf²とすると、conic関係から
`q=(R-A²)/(S-A²)=1+(K-1)R/A²`。
tauはu→-u,K→1/K,R↔Sなのでq→1/q。
従って `w=u(q-1)/(q+1)` はtau固定で、
`R=2A²w/((K-1)(u-w))` が逆式になる。
F5=E5(u)、u²∈E5を用いると `F5(R)=E5(w)(u)`、tau固定体はE5(w)。
tauはV4の二つの独立符号を交換するため、この有限source作用は合計order8となる。

Sp4のactual frame hを加えた元source fieldをS、actual full fieldをL=k(Y16)とする。
上のfull deck同値から `L⊂k(Sp4)E5(w)`。
右の固定体はS内index8である。
既受理FULL_SOURCE_DEGREEは、全normal form代表4つ×Aut component2つを数え、
**actual full source map**のdegree8を証明している。
普遍ordered spectral曲線のm-degree8を流用したものではない。
両indexが同じなので包含は等号。既受理の同じbase embedding
`k(B15)=k(Sp4)E5` により `L=k(B15)(w)` を得る。
QQ上の有限有理同値とchar0 degreeを使い、CからQQへ逆降下していない。

## Jのcompact式と次数

独立に一次式L(w)=L0+L1wへ整理すると

```
L0=u(K-1)(-Iv*A²*K+8u(K+1)),
L1=Iv*A²*K(K+1)-8u(K-1)²,
J(w)=L(w)²/[K*u²*(K-1)²*(u²-w²)]。
```

保存されたJ(w)全式と、元vを回収したJ(t)の両方に一致させる。
L(±u)は一般に非零なので分母との相殺は無く、actual degreeは2。
w(t)のdegree4との合成でJ(t)のdegreeは8となり、実際の約分後次数も確認する。
tau作用Iv→-K Ivも五定数から独立に確認し、上式がtau固定であることを照合する。
係数を旧sectionの名前の同じ別函数と同一視しない。

full150再展開、点、prime、rank、不変式source26、author rerunは無い。
二重被覆の新v0座標表示は別の同時review receiptへ保存する。
