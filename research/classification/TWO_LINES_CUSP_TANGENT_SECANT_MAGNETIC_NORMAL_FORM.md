# 尖点接線・割線源の全磁気正規形と一つの全parameter認識

2026-11-07、PRIMARY `E_reduced_nonnormal`。全六磁気源の可逆な還元を先に固定する。
この還元だけでは既知成分への所属は主張しない。

## 全源と磁気translationの逆

baseをh,x,y、fiberをz,wとする。対象は
`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json` の
`cusp_cusp_tangent_cusp_secant` で、C=Y²Z−X³、二直線Y,Y−Xである。
全括弧は

```
{x,z}=−hx,  {x,w}=3hx−3x²,
{y,z}=3(h−x)²−3hy,  {y,w}=−9xy,
{z,w}=3xz+3(r0h²+r1hx+r2hy+r3x²+r4xy+r5y²),
{h,−}={x,y}=0.
```

U=−hx∂x+[3(h−x)²−3hy]∂y、V=(3hx−3x²)∂x−9xy∂yと置く。
z'=z+ell、w'=w+kappaという実PGL translationは、磁気多項式を
`R→R−3x ell−U(kappa)+V(ell)` へ送る。
ここで

```
ell=(r3+r2)x/2+r4y/4,
kappa=[−3r1−3r3/2+9r2/2]x−r2y
```

を取ると、新しい磁気項は正確に `3(Jh²+Iy²)` となる。
I=r5、J=r0+r2であり、全六係数は一つも反転していない。
逆変換はz=z'−ell、w=w'−kappaで、baseは不変だから同じ多項式式である。
したがって全sourceに正則な座標変更であり、特殊な磁気点だけの正規化ではない。

I≠0ではfiberを `(z'',w'')=(z'/I,w'/I)` とし、pushforwardしたtensorをI倍する。
全anchorとfiberの一次項は元に戻り、磁気項だけが

```
R_lambda=3(lambda h²+y²), lambda=J/I
```

となる。逆fiber変更と逆tensor scaleも可逆である。
従ってI≠0の全六源のPGL/scalar帰属は、この一つの全lambda源の帰属と同値である。

## I=0の閉補集合を戻せる条件

元のA6 sourceは既約で、I≠0は稠密開である。
既知成分の有限閉和Kに全lambda源が入ることを別に証明できれば、
PGL/scalar不変性で元のI≠0全像もKに入り、その閉性でI=0も全てKに入る。
これはI=0を反転した同型をそのまま代入する議論ではない。
逆に一般lambdaがKの外なら、I=0の帰属又は全特殊lambdaの帰属は別に残す。

## 全lambdaに同時に有効なDarboux chart

rank4評価点を `(h,x,y,z,w)=(1,1,1,0,0)` に固定する。
全lambdaで実Pfaffianは−9である。x,yをそれぞれ1だけ平行移動した後の評価行列に対し、
四行のDarboux frameは

```
e_x, −e_z, e_y−3e_x, ((lambda+1)/3)e_x−e_w/9
```

である。detは−1/9、全lambdaについて可逆で、source parameter分母は無い。
したがって一つのQQ[lambda]計算でfaithful S61の全150逆再構成、全134式、元の全350式、
T39/T40の各全45graph式を検査できる。

両graphの全残式はQQ[lambda]のidealとして保持する。
全て恒等零なら全lambda所属、そうでなければ全残式のgcdの零点が正確な特殊所属集合である。
定数gcd1なら例外lambdaも無い。parameterの無限遠はI=0側として元sourceに保持する。
global成分とS61全graphの同値には既受理
`S61_GLOBAL_COMPONENT_EVALUATION_FIBER.md` とその独立reviewを用いる。
陰性を他の既知親又は全ambient成分の不存在へ読み替えない。

## 有限実行の停止規則

既点の全六法線二次mapは非孤立と確定済みで、再実行しない。
今回の計算は全六磁気源を保つ一つの可逆正規形と全lambda所属判定だけである。
一つの全graph計算の後にanchor・frame・高次jetを追加しない。
