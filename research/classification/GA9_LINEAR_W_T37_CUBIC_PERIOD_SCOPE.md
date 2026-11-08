# Ga9 linear-w と T37：全 cubic-period の恒等式と局所 flag の射程

2026-11-08、PRIMARY E。凍結済み true completed T37 と既受理の全源定理を手で読解した純数学記録。新しい point、operator、normal/profile、jet 延長又は CAS は実行しない。Ga9 の T37 帰属・非帰属、新成分性は未解決のまま。

## 1. T37 の一般軌道には order二零点がない

既別作者 `K4A_ORDER_TWO_ZERO_SEPARATE_REVIEW.md` は、K4a actual24 親の全 intrinsic 条件 `∂(E∧P)(q)=0` を50×15行列として再構成し、既指定十五行 minor の determinant

```
3411634563591564165120 ≠ 0
```

を得た。従ってこの同じ dense24 軌道には order≥2 の projective zero が存在しない。後の `K4A_T37_ACTUAL_PGL_IDENTIFICATION.md` と別作者全150係数恒等式は、この軌道を true completed T37 と同定している。よって T36 の distinguished order二零点における quadratic-pencil の閉条件を T37 へ移すことはできない。Ga9 にはその零線上で全 firstjet が零であるが、それだけでは T37 の全閉包との比較にならない。

## 2. 全 canonical T37 pair に残る product identity

`CANONICAL_INDEX_TWO_NODE_WHOLE_ORBIT_ENTRANCE.md` の全可逆正規化後、

```
q=ac−b², A=c²+bd, B=a²+be,
C=(AB−q²)/b=c²e+a²d+2acb+bde−b³,
F=det[[C,A,B],[A,a,b],[B,b,c]],
g1=bB−cA, g2=bA−aB.
```

除算記号を含む C は右辺の literal polynomial として使い、b=0 の空間 chart を捨てない。全五座標で恒等的に

```
g1*g2=q³+b*F.                                           (1)
```

実際 F=Cq−cA²+2bAB−aB² だから、

```
bF=q(AB−q²)−bcA²+2b²AB−abB²
   =(ac+b²)AB−bcA²−abB²−q³
   =(bB−cA)(bA−aB)−q³.
```

全積の恒等式であり、F=0 上、b unit chart 又は正規化の一般点だけの関係ではない。従って任意の actual moving linear frame とその逆にも保持される。source の residue は λ*dlog(g1/g2) だが、二つの cubic と q、F を自由に別々に Taylor 化することはできない。既受理 fixed-D 認識へつなぐ actual family は、この関係を保ち、全 tensor の負冪を消す必要がある。

## 3. Ga9 の H-component residue は cubic confluence と整合する

canonical 四倍 Ga9 を h=0、y=1 chart の (X,Z,W)=(x/y,z/y,w/y) へ制限すると

```
{X,Z}=3Z+4W², {X,W}=−3W, {Z,W}=4.
```

その三次元 volume に対する contraction は、非零定数 convention を除いて

```
dC, C=4X+3ZW+4W³/3
     =(12xy²+9zwy+4w³)/(3y³).
```

三括弧で C の Casimir 性を直接確認できる。h*F のもう一方の因子はこの chart で h=0 上1となるので、inverse の H-residue も同じ dC である。

N=(12xy²+9zwy+4w³)/3 とすれば、形式的な cubic-pencil identity

```
ε^(-1)*dlog((y³+εN)/y³) → d(N/y³)
```

は成立する。したがって「Ga9 に三重分母の exact residue がある」という degree 情報だけで T37 を排除できない。しかしこの単独 identity は (1) や ambient Poisson の正則性を含まないので、actual family、帰属又は endpoint tensor を証明するものでもない。

## 4. secant endpoint の global contraction も pencil にはならない

true completed T37 の全 period を1とし、座標を (H,x,y,z,w) とする。secant endpoint q=[0:1:0:0:0]、x=1 chart の intrinsic firstjet は

```
H*∂y∧∂w.
```

非零なので quadratic part 単独は intrinsic leading tensor ではない。さらに全 cubic `E∧P` の volume dual 二形式を q で contract して得る一形式 α についても、係数は

```
α_H=2yzw−H²y/2−Hxz−z³/3+Hw²/3,
α_y=2Hzw+H³/2,
α_z=−3Hyw+H²x+Hz²/3+w³,
α_w=−Hyz−H²w/3−zw², α_x=0.
```

全体の符号は volume convention にのみ依る。x 依存部は `xH(H dz−z dH)`。従って α∧dα の dx∧dy∧dz 係数は

```
−H³*(2zw+H²/2) ≠ 0.
```

この distinguished endpoint の global contraction 自体が非可積分である。局所 firstjet を捨てた quadratic part 又は global contraction を、T36 と同じ rational-pencil foliation とみなすことはできない。

## 5. 次の whole 実族に必要な条件

有効な positive criterion は既 fixed-D 全源定理のまま：真 T37 に属する actual curve が regular な非零 endpoint を持ち、その signed divisor が Ga9 の hF に一致し、unique singular plane 上の rank が二なら、全磁気逆写像により Ga9 unit 全23軌道を認識できる。(1) はその全 cubic-period 源を保つ必要恒等式である。divisor 又は residue の一致だけを actual tensor の一致へ置き換えない。

本記録は order二 flag の適用不能と全 product identity を明示する。whole T37 negative、Ga9 の新成分又は既知親への帰属を追加しない。全九 SL2 ansatz の一回陰性は別の凍結記録のままで、補正項を延長しない。global 台帳、主稿、code/current contracts、Git は変更しない。
