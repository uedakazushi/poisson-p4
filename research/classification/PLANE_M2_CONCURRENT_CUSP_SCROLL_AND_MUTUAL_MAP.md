# concurrent double残二軌道のcusp-scrollと実相互同定

2026-11-08、PRIMARY E、純座標代数。既fullA8／全normalizerの後、同じcanonical二点だけを読む。ambient profile／GB／新しい点を試さない。以下のliteral証明書は後続jobで別に保存する。

Aはperiod5/2、Bは−1/2、各canonicalは保存fullsourceのp1=1、他P1係数零であり、P2は全零。実五次式F_A,F_Bをsigned Euler-Pfaffianとして取る。

## 1. Homogeneous cusp恒等式

```
Q_A=9au−2v²,
C_A=16cu²+9(a−4b)uv−2v³,
486 u F_A=−2 C_A²−Q_A³.

Q_B=9au+2v²,
C_B=48cu²+(27a−36b)uv−2v³,
162 u F_B=2 C_B²−Q_B³.
```

u一chartでは(a,c)を次のregular coordinatesへ置換できる。

```
A: F_A=−256/243[c+9(a−4b)v/16−v³/8]²
       −3/2[a−2v²/9]³,
B: F_B=256/9[c+9(a−4b/3)v/16−v³/24]²
       −9/2[a+2v²/9]³.
```

従ってこのchartの全singular supportはcuspの二座標零であり、smooth cubic scrollのdense開である。その全閉包S′は、次の2×3 matrixの全minorsで定まる。

```
A: [[u,2v,9b],[v,9a,4c]],
B: [[u,2v,3(a−b)],[v,−9a,−4c]].
```

row二を二倍して変数を並べると標準行列[[y0,y1,y3],[y1,y2,y4]]となる。これはS(2,1)のsmooth scrollである。u零ではFのa偏微分が非零定数×v四となるので、残るsingular点はu=v=0だけ。全偏微分はそのplane S上で零である。従って**reduced SingFはS∪S′**であり、scrollだけとしない。scroll境界の二linesはS内に含まれる。

uchartのintrinsic bracket（raw一成分ではない）を使うと、S′上で

```
A: {b,v}_affine=4(9b−v²)/27,
B: {b,v}_affine=4(9b+v²)/9.
```

従って両scroll上のrestriction rankは一般二、generic transverse特異点はA2 cuspである。元plane Sのrank二restrictionの零cubicはab(a−b)、三直線はconcurrentのままである。

両Fはcについてprimitive quadraticで、discriminantのv零restrictionはそれぞれ−512a³u⁵/81、512a³u⁵。これらが非平方なので絶対既約である。五偏微分のconstant relationでは、c²u²,c u³,a²u²,b²uv,buv²の係数を順に比較すると五定数が全て零となる。従って本質変数は五であり、直接essential三定理へは入らない。

## 2. 二canonical tensorは実際に同じ軌道

```
M_AB:(a,b,c,u,v) -> (−a,b−a,c/3,u,v), det=−1/3.
```

逆は(A,B,C,U,V)→(−A,B−A,3C,U,V)である。十括弧の直接係数比較は

```
(M_AB)_* Π_A=(1/3)Π_B
```

というRAW恒等式を与える。Euler補正は不要である。例えばuv、au、buの三raw括弧はそれぞれ同じ1/3比となり、base及びc-rowも同じ線形置換で一致する。後続literal jobでは十括弧全150係数を全て検査する。

Q_B(M_AB y)=−Q_A(y)、C_B(M_AB y)=C_A(y)なので

```
F_B(M_AB y)=−3 F_A(y),
D_{M_*Π_A}=(1/9)F_B.
```

従って二つの全A8族は同一のcanonical実軌道閉包へ入る。二つの独立なambient profileを追加する理由は無い。現段階でその軌道閉包の全ambient成分性又は親の名前は未判定である。

## 3. 同一smooth-scroll divisorへの標準frame

repo標準scroll quadricsを

```
q1=x0x3−x1x2, q2=x0x4−x1x3, q3=x2x4−x3²
```

とする。次の実framesはscroll S′をこの標準scrollへ送る。

```
A: X=(72b,64c,u,2v,18a),
B: X=(72(b−a),192c,u,2v,−18a).
```

BのframeとM_ABの合成はAのframeと同じである。両Fは同じ固定五次式

```
F*=x2 q1²−2x3 q1q3+x4 q3²
```

に移り、正確に

```
3888 F_A=−F*(X_A),    1296 F_B=F*(X_B),
x2 F*=(x2 q1−x3 q3)²+q3³.
```

この固定divisorは旧rank二smooth-scroll A25 cuspのprimitive first-row型と一致する候補である。第二row／closed inverse residueとglobal Poisson extensionを戻さず、divisorの同一性だけからtensorの既親T36–T40を認定しない。次のrecognitionはこの実frame／全fixed-D源の構造を用いる。
