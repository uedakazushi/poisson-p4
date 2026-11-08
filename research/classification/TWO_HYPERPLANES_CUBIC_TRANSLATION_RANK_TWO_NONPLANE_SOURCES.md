# rank二平行移動の二超平面・三次型：非平面支持の二つの完全源

2026-11-08、PRIMARY `E_reduced_nonnormal`。保存済み全Jacobi係数の手による消去。
二つの全源を明示するが、周囲成分性・既知親・新成分数はまだ主張しない。
既点のjetやnormal mapは再実行しない。

## 1. 全係数入口と保持する平面補集合

`TWO_HYPERPLANES_CUBIC_SUPPORT.md` §5.1の最終型
`V=x∂u+h∂v` を扱う。h-Casimir代表の全括弧は

```
{x,u}=x², {x,v}=hx, {x,w}=0,
{u,v}=h²n0+hn1x+hn2w+(1−n4)hu+n3x²+n4vx+n5wx+n6w²,
{u,w}=h²n7−hn11u+hn8x+hn9w+n10x²+n11vx+n12wx+n13w²,
{v,w}=h²n14+hn15x+hn16w−hn18u+n17x²+n18vx+n19wx+n20w².
```

これが全21自由係数の線形源であり、残るJacobiの全係数は
`CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json` の
`j111_rank2.remaining_coefficient_equations` に固定される。
三次因子はF=−h P_uw+x P_vwで、実五次はhxF。

T=hu−xv、A=n11、B=n18と書くと、Fのu,v部分は `(Ah−Bx)T`。
A=B=0はu,v非依存の本質三変数枝であり、既存の平面十三配置側にそのまま残す。
以下は(A,B)≠(0,0)だけを尽くす。先に全sourceの低次係数を零に置かない。

## 2. 全可逆normalizerと二つのfinite chart

A≠0ならw'=w/AによりA'=1、B'=B/Aとなる。
Vはwを含まないので、この変換はVをそのまま保つ。
r=B/A≠0ならさらにx'=r x,u'=r uを施すとA'=B'=1で、
`{x',u'}=x'²,{x',v}=hx'` も保たれる。
従ってA≠0は(A,B)=(1,0)又は(1,1)に尽くされる。

A=0,B≠0は二超平面を交換して同じchartへ移る。
具体的に `(h',x',u',v',w')=(x,h,v,u,w)` とすると、pushforwardしたH行は
h'で割れ、Euler gaugeでh'-Casimir代表へ戻してtensorを−1倍すれば同じV正常形となる。
三次の非平面項は

```
(Ah−Bx)(hu−xv) → (B h'−A x')(h'u'−x'v'),
```

なのでA,Bは交換される。これはprojective bivectorの可逆PGLと代表変更だけである。
A=B=0は一度も反転せず、別の平面補集合として保持する。

## 3. Jacobi消去の共通骨格

A=1,B=rとする。保存された全係数式から順に

```
n16=r n9, n19=r(n12+1), n20=r n13,
n6=n13(n4−1),
n2=−2n13n14+2r n13n7+(n4−1)n9,
n0=−n14n9+r n7n9+(n4−1)n7,
n1=r n0−(n12+n4)n14−n15n9+r n9n8+r(n12+1)n7+n4n8,
n5=n12n4−2n13n15−r n9n4+r n2,
n3=r n1+r n10n9+(n4+1)n10−(n12+n4+1)n15−n17n9+r(n12+1)n8
```

を得る。最後の係数式は厳密に−2r n13となる。
残る二係数式と合わせ、r=0,1の全解は次の表示になる。
導出にr以外の未知係数の反転は無い。r=0で必要なn17非零だけは既約三次の条件から正当化する。

## 4. L=h：完全八parameter源

自由変数は

```
n4,n7,n8,n9,n10,n14,n15,n17, n17≠0.
```

残る全係数を

```
n11=1, n18=0, n12=−n4−2,
n13=n16=n19=n20=n6=0,
n0=−n14n9+(n4−1)n7,
n1=2n14−n15n9+n4n8,
n2=(n4−1)n9,
n3=(n4+1)n10+n15−n17n9,
n5=−n4(n4+2)
```

とする。逆にこの表示は保存された全Jacobi係数を恒等的に零にする。
必要性の最後の二式は
`n17(n12+n4+2)=0, n13n17=0` である。
実際F|h=0=n17 x³なので、n17=0ならF全体がhで割れる。
対象の既約三次を保つためn17≠0が必要であり、その補集合を黙って捨てていない。

全三次因子は、可逆linear変更

```
U=u−n9w−n7h+(n14−n8)x,
V=v−(n4+2)w−(n15−n10)x
```

により正確に

```
F=h²U−hxV+n17x³
```

となる。これらの変更はcofactor V=x∂u+h∂vを保つ。
さらにh,x,u,v,wの対角尺度とn17の非零平方根によりx³係数を1へ送れる。
全域でprimitiveなU一次式だからFは絶対既約で、本質変数数は正確に4である。

## 5. L=h−x：完全七parameter源

自由変数は

```
n4,n7,n8,n9,n10,n14,n15
```

の全A7。残る係数は

```
n11=n18=1, n13=n20=n6=0,
n12=−n4−1−n9, n16=n9, n19=−n4−n9,
n17=n10+n14−n7,
n0=−n14n9+n7n9+(n4−1)n7,
n1=n14−n7−n15n9+(n4+n9)n8,
n2=(n4−1)n9,
n3=n14−n7+(n4+1+n9)n10−n17n9,
n5=−(n4+1)(n4+n9).
```

これも全Jacobi式と同値で、自由parameterの零点を一つも除外しない。
L=h−xとして、全Fは

```
F=L(hu−xv)+[hx+L(−n9h+(n4+n9)x)]w+C_base(h,x)
```

である。u'=u−n9w、v'=v−(n4+n9)wとすればwとの余分な項は吸収される。
`C_base=c0h³+c1h²x+c2hx²+c3x³` と書き、さらに

```
U=u'+c0h, V=v'+c3x,
W=w+(c1+c0)h+(c2+c3)x
```

と置くと、厳密に

```
F=(h−x)(hU−xV)+hxW.
```

これは全parameterに正則でdet1のlinear変更であり、cofactorも保つ。
U,V,Wに関する三係数h(h−x),−x(h−x),hxのgcdは1なので、primitive linear polynomialとして
全域で絶対既約である。三係数はSym2(h,x)を張るため、本質変数数は正確に5。
この三次はPerazzo cubicのlinear正常形であるが、その外部分類定理は使用しない。

## 6. 未完了義務と有限な再開対象

ここまででrank二平行移動の非平面全支持は二つの明示的な完全源と全可逆normalizerへ尽くされる。
平面A=B=0、可約三次、rank低下点はそれぞれ明示的に別義務又は対象外として保持した。
全parameterで元の21係数sourceを用いるので、任意の低次項を先に消した部分sliceではない。

二源の周囲成分性・既知親は未解決である。
新成分数を増やすには、真のambient全D、実全源下界、必要最小の支持証明、独立reviewと全prior比較が別に必要。
既知MPS等へ入れるには、全源を保つactual curve又は別の全支持定理が必要である。
本稿の有限正常形及び三次因子のgeometryだけからその結論へ飛ばない。
