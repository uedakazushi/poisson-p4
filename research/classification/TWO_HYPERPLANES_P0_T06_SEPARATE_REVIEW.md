# j21(h+x,h)・p零の全五源：named T06 confluenceの独立確認

2026-11-08、PRIMARY E、`/root/log_manuscript`。**全五源と全PGL閉包のT06帰属PASS**。
producer `TWO_HYPERPLANES_P0_T06_CONFLUENCE.md` の最小literal曲線を別作者として確認した。
既独立Fraction checkerのisolated `--p0-t06` 一回で、実測0.147762959秒、RSS26940 KiB。
producer、S61 engine又はCASをimportせず、新point、anchor、normal、jet又はGBを使っていない。
68未知数matching、cokernel rank又は探索profileの完全性は被覆の依存にしない。

## 1. 全五sourceのnormalizerは全parameterで可逆

`TWO_HYPERPLANES_CUBIC_SUPPORT.md` §5.2とaccepted dataの
`j21_hx_h_p0` を必要十分な全source入口として用いる。
元十四自由係数源へaccepted bindingを戻し、別に表示五sourceからも組み立てて全150一致を確認した。
全strict10、元350、旧全remaining equationsも恒等零である。

free parametersはn2,n3,n5,n8,n10。
A=n10*n3+n5+n10²/2、B=n2+n3*n8として

```
U=u+Ah-Bx, V=v, W=w-n10h-n8x
```

を独立に組み立てた。det一、全inverse、両向き全150係数を確認した。
source parameterの反転は無い。l=n3とすると全sourceは

```
xu=xv, xv=x(h+x), xw=hx,
uv=hu-v²/2+l[-hv+(h+x)w],
uw=hu-w²/2+lh(w-v), vw=h(v-w)-lhx,
D=hx F0, F0=[h(v-w)²+xw²-2hxu]/2
```

となる。F0はuについてprimitive linearで絶対既約、全lに依存しない。
四つの余分な係数をchosen pointへ置いたものではなく、全逆付きwhole正規化である。

## 2. 真のrawT06全四parameter源

accepted raw sourceは座標x0,...,x4、free p0,p1,t0,t1について

```
P01=t0*x4², P02=t1*x3², P12=p1*x1*x2,
P13=-p0*x1*x3, P14=(2p0+p1)*x1*x4,
P23=-(2p0+p1)*x2*x3, P24=p0*x2*x4,
P34=(3p0+p1)*x3*x4/2
```

で、他二entryは零である。
この全四parameter表示を独立に組み立て、rawT06 packetの全十括弧と一致させ、元350式を恒等零と確認した。
このaccepted actual sourceのgeneric period/edge開はT06に入る。
全四parameter写像はpolynomial、T06は閉なので、同じ全源のparameter境界もT06に入る。
新しいSNC一点だけから任意の別族へ帰属を逆向きに推論しない。

## 3. 四literal profileと全PGL inverse

親の真のparameterを

```
p0=e+2e²l, p1=-e-4e²l, t0=e, t1=-e
```

へ同時代入し、e非零でactual PGL

```
x0=-2u-2v/e-(h+x)/e²,
x1=x, x2=h, x3=w+h/e, x4=v-w+x/e
```

を行う。det負二と全逆を独立確認し、余分なtensor scaleを付けない。
source parameterの積を独立なTaylor係数へ緩めていない。

全quadratic liftのcanonical Euler代表を

```
C(P)=P-E wedge div(P)/5
```

とする。div(E wedge L)=5L-(trL)Eなので同じprojective tensorを表す。
CのkernelはEuler vertical directionsそのものであり、都合の悪いPole行を捨てる操作ではない。
全150にこのprojectionを施し、divergence零とactualD不変を直接確認した。
rawT06全parameterのcanonical代表がstrict Poissonとなることも別に全十式で確認した。

literal曲線のC(P)全150係数ではe負冪が消え、QQ[e,l]の多項式となる。
全e,lでstrict10及び元350が零。
e零の全150は§1模型のC(P_l)へ厳密に一致し、同じprojective tensorである。
normalizerのinverseにより元全五sourceへの端点一致も全parameterで戻る。

Euler-signed actualDを両方のtensorから直接作ると

```
D(e,l)=(1+el)(h+ew)(x+e(v-w))*F0
```

に完全一致し、親五次のactual pullbackにinverse determinant負1/2を掛けた式とも一致する。
単なるdivisor一致でPoisson帰属を推測したものではない。
lを反転せず零を含む全lが同じ正則端点を持つ。
e非零は真のT06源のinvertible PGL像、e零は正則極限だから全endpointがT06内にある。
T06のPGL安定性と閉性、全source normalizer inverseから、全五源及びその飽和閉包の帰属が従う。

## 4. 保存とscope

reportは `CLASSIFICATION_TWO_HYPERPLANES_P0_T06_SEPARATE_REVIEW_20261108.json`、
jobは `two_hyperplanes_p0_t06_separate_review.json`、
入口は `inputs/two_hyperplanes_p0_t06_separate_review_manifest_20261108.json`。
旧checker全FunctionDef ASTと五currentjob exact bytesは
`inputs/two_quadrics_snc_pre_p0_t06_review_source_contracts_20261108.json`へ保存した。
全producer、旧source/report/current契約は不変である。

新成分、局所一意性、ambient smoothness、他nilpotent型又は全Eの被覆は主張しない。
残三plane1+1+3 source、Ga9linear-w、stable二源などは別の未完了義務である。
