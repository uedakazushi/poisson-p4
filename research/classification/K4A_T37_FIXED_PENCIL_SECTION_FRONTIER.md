# T37 固定 P1 の全九 equivariant plane：全sectionと正確な停止範囲

2026-11-08、PRIMARY E。`K4A_T37_FIXED_PENCIL_COISOTROPIC_PLANES.md` の純手分類を使う。
新有限五 flag の一括 exact run は **PASS、6.136秒**。
全 actual coordinate 行列と逆、十括弧、coisotropy、anchor、pole-free normalcone、
端点の全350 projective Jacobi、全六列の磁気 defect を
`reports/CLASSIFICATION_K4A_FIXED_PENCIL_COISOTROPIC_PLANES_20261108.json` に保存した。
共有 runner を変更していない。独立再算術は未実施であり、新componentを計数しない。

## 有限五 flag の全section

base 座標を `(C,B,V)` とする。非零定数を除き、I± の section は `Q2*G3`。
r²=6、符号を r の二値に含めると

```
34Q=34B²+158rBC-40rBV+1644C²-690CV+75V²,
G=B³-r*B²V/2+72BC²-6BCV+3rC³-3rC²V.
```

Q の対称行列の determinant は -67500 で smooth irreducible conic。
G の V 係数は `-r*(B+rC)²/2`、V 非依存項の B=-rC への制限は `-75rC³`。
従って両係数は C[B,C] で互いに素、primitive な V 一次式 G は C 上絶対既約。
この全二 flag は [2,3] 型であり、二重 H 又は二つの simple H を持たない。

II± は r²=2 として

```
T=B+33rC/5-9rV/10,
X=B+3rC,
D=const*T²*[T*X²+(81/5)*C²*(2X-9rC)].
```

括弧内は T 係数 X² と非依存項が互いに素なので絶対既約 cubic。
唯一の cusp は X=C=0。affine T=1 の初項は `X²-(729r/5)C³` で A2。
H=T=0 は cusp を避け、交点 `C²*(2X-9rC)` は ordinary tangency と一単純点。
これは既受理 nil22 image double_A 配置で、残 stable 二配置ではない。

infinity flag `b=u=0` の section は

```
T=8B+3C-3V,
D=const*T²*[T*B²+9*(B-2C)²*(B-3C)].
```

同じ primitive T 一次式判定で cubic は絶対既約。
B=C=0 の唯一の cusp を T=0 は避け、B=2C で二重、B=3C で単純に交わる。
従ってこの flag も double_A 配置。

## 全九 flag と全磁気補集合

旧 k=0 四 flag の全section は、別の不変証明書から

```
C²*(V-9C)³;
(3V-19C)*[B²+2C*(3V-19C)]²;
C²*(cuspidal cubic), 二符号、H は cusp を避け ordinary tangent+secant。
```

新四＋infinity一と合わせ、全 U=0 固定成分からの指定正則 equivariant plane class は
正確に九 flag で尽くされる。どの section も
`H²*(line*smooth conic)`、`H²*(three distinct lines)` 又は
`H1*H2*(irreducible cubic)` の残配置を与えない。

新五 flag で全 first-order plane motion 六列の磁気 defect はいずれも rank 五。
全6×6行列、rank五 minor、全 pole-free normalcone tensor を保存し、
一つの residual magnetic quotient を捨てていない。
divisor 一致だけから full sixmag source membership を推論しない。
この場合、section の因子型が既に目的と異なるので residual quotient の延長は行わない。

## 陰性の適用範囲

これは固定 P1 の全正則 equivariant plane class と uniform two-normal contraction の結論。
正則 base frame と arbitrary first-order plane motion を加えても section は同じなので、
この class だけでは残 stable 二源・残 plane1+1+3 三源を覆えない。
全 T37 の非所属、別固定成分、singular base frame、nonuniform contraction 又は
真の moving-source confluence の非存在は主張しない。
同じ class の profile や点を増やさず、別の whole-source geometry へ切り替える。
