# 普遍平面四次の二節点とgenus上界：独立レビュー

2026-11-04。**二つの通常節点と、各幾何的商成分の正規化genus≤1を受理する。**
対象はordered-eigenvalue baseと既存relative-sign quotientの後の平面四次Psi。
元の二次coverまたはfull残余曲線のgenus、exact genus1、全四次の既約性は主張しない。
authorの別の(2,2)模型のsmoothness計算は今回の監査に含めない。

## 独立に再構築した係数恒等式

checkerはauthorを再実行せず、独自の小さいsparse Fraction多項式積で

```
W=2uvD+(u+v)rt, m=u-v,
4uv*Psi=W²-m²(r²+4fuv)(t²+4cuv)
```

を全係数照合する。m=0、u=v=u0でD=B00+2ar+2dtとすれば
`Psi=P0²`、`W=2u0 P0`、`P0=rt+u0(B00+2ar+2dt)` も独立に照合する。
凍結済みactual quartic tableの全次数4、r⁴係数16x⁴y⁸、15係数/1592項のbindingを確認する。
元の1592項の新規展開やM2は行わない。

係数体は `K=QQ(x,y,a,k1,k2)`、`K0=QQ(x,y,a,k2)`。
`K10=k1/2-A0|m0` はK0上超越的であり、u0,B00,d,c0,f0はK0に属する。
xy≠0かつu0≠0である。後者はu0の分子のk2係数が1であることから従う。
既受理のgeneric spectral chartを使い、ある一つの係数点での性質を一般へ外挿しない。

`t=(K10-xr)/y` を代入したP0のr²係数は−x/y。
yを払った多項式は

```
yP0=-xr²+(K10+2u0ay-2u0dx)r+u0B00y+2u0dK10。
```

そのdiscriminantのK10²係数は1、P0自体では1/y²で非零である。
従って幾何的係数体で相異なる二根を持つ。QQ上に二根が有理的にあるとは仮定しない。

## 各根で二つのqが単元である理由

qR0=r²+4f0u0²がP0と共通根rを持つと仮定する。
qR0はK0上のmonic polynomialなのでrはK0上代数的。
K10はその有限拡大上でも超越的だから、P0のK10係数を零とする必要があり、
`r+2u0d=0` となる。すると
`qR0=4u0²(d²+f0)`。独立なdenominator-cleared照合により

```
d²+f0=2ax+3x²+y² !=0
```

なので矛盾する。この議論でrがK10を含む一般の根であることを勝手に否定していない。
共通根という追加仮定こそがK0上代数的にする。

同様にqT0=t²+4c0u0²を仮定するとtはK0上代数的となり、
`r=(K10-yt)/x` をP0に入れた式のK10係数から `t+2u0a=0`。
しかし `a²+c0=-2ax-x²+y² !=0` なので矛盾する。
よって二根の各点でqR0*qT0は非零の単元である。

P0の二根は単根なので `∂W/∂r=2u0 P0'(r)` は非零。
従って(m,W)はplaneの局所正則座標で、曲線の二次初項は
`W²-qR0*qT0*m²`。標数0の幾何的係数体上で相異なる二直線に分かれ、
各点は通常節点、delta=1である。他の特異点が無いとは言わない。

## genusと失敗したprobeの範囲

幾何的に既約なdegree4 supportなら算術genus3から二つのdeltaを引き、正規化genus≤1。
可約またはnonreducedなら、各既約supportのdegreeは高々3なので同じ上界を持つ。
従って全ての幾何的**平面商成分**に上界1が成立する。
これをh²=qRで与えられる元の二次coverのgenusへそのまま移してはいけない。

旧係数点(1,2,3,5,7)/32003でのprojective smoothness testは
`all_passed=false` のままhashで保つ。新しい点、prime、旧テストの再実行は無い。
失敗した特殊fiberからgeneric singularityを推論していない。
今回の二節点は普遍恒等式と超越係数による別証明である。

再現:

```
python code/q_firstjet_residual_forced_nodes_review_20261104_v63.py \
  --output /tmp/universal-quartic-two-node-review
```
