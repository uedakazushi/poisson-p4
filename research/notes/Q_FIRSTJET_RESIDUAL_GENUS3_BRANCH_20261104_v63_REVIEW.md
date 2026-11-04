# Proper二節点模型とfull ordered曲線の種数3：独立監査

2026-11-04。**proper商曲線の幾何的種数1、full ordered二次coverの幾何的種数3、
実際の選択成分のm-degree8を受理する。** 対象は普遍W16像のspectral fiberについて、
明示したordered-eigenvalue baseと既存relative-sign quotientの後の曲線である。
NEW hypersurfaceのmultisection degreeやWeyl quotientのgenusへ広げない。

## 同じ係数点のproper模型を独立にbindする

既受理の二節点証明と、元のactual1592係数・全algebra逆写像を監査した
RECIPROCAL_GENUS_REVIEWを再利用する。元のsmoothness失敗はfalse reportのまま保つ。
今回の専用checkerはauthorを再実行せず、新しい点・prime・GBも使わない。

同じ(1,2,3,5,7)/32003の15係数を、独立な2×2 companion-matrix算術で
`2alpha²-287alpha+10914=0` 上へ移す。beta=287/2-alpha、delta=alpha-beta。
保存されたpair乗算とは別に行列を掛けて、

```
delta^4 F((alpha V-beta U)/delta,M,(V-U)/delta)
```

のM=1での全九係数を再構築した。U,V次数が各2以下、deltaとU²V²係数のnormが
非零であることも照合する。四つのpatchでの全係数を保存M2 scriptへ個別にbindし、
四unit log、正常exit、全hashを確認する。M2を再実行したとは主張しない。
alphaの一根だけを選んでおらず、scriptは二次étale代数の全幾何的根を扱っている。

## Proper/flatな移行とgenus1

Z_(32003)上の係数baseを、x,y、節点discriminant、必要な係数分母で局所化する。
monic P0の根を付けた有限étale base上でdeltaは単元。
`F(r,0,z)=lambda(r-alpha*z)^2(r-beta*z)^2` なので、変換後のU²V²係数は
`lambda*delta^4`、lambda=16x⁴y⁸で単元となる。
単元係数をleading monomialにしたmonic divisionで各bigraded pieceはbase-free。
従ってこの(2,2) divisorはproper flat familyである。

相対non-smooth locusのproper imageは閉。特殊baseの二根全てのfiberで四patchがsmoothなので、
baseをさらに縮めると全node choicesでsmoothになる。この非空openはchar0 generic点を含む。
有限étale coverがsplitしている場合も、全rootを検査し有限proper imageを除くことで同じ結論。
一つの非proper affine fiberからのemptinessやnilpotenceの外挿ではない。

幾何的smooth (2,2) divisorはconnectedである。例えば
`H^1(P1×P1,O(-2,-2))=0` とdivisor完全列からH0(O_curve)=kが従う。
smooth connected曲線は既約で、adjunctionによりgenus=(2-1)(2-1)=1。
平面の二節点をblowupし、結ぶ直線をcontractする際、strict transformとの交点は
4-2-2=0、直線は曲線成分でないので曲線から離れている。
従ってproper (2,2)模型が元四次の正規化を与える。
元の幾何的既約四次の算術genus3との差は既知二節点のdelta計2で尽くされ、他に特異点はない。

## Finite部での分岐を全て除く

K=QQ(x,y,a,k1,k2)、K0=QQ(x,y,a,k2)とする。
qR=r²+4fuv、qT=t²+4cuv、t=(K1-xr)/y、K1=k1/2-A0(m)。
全係数と同じbaseを保持し、k1の超越性を使う。

独立sparse Fraction算術で、basis(1,r,t,rt)、r²=-4fuv、t²=-4cuvにおける
Wの全multiplication matrixとdeterminantを再構築する。
二回のnormを**逆順に**取った式とも比較し、authorの

```
Norm(W)=16(uv)^4(E²+64fuv F²),
E=B0²-16fuv*a²+16cuv*d²-4cf(u+v)²,
F=aB0+2cd(u+v)
```

へ一致させる。実際のm-leading dataからEのdegree4係数は
`4-4(-(x²+y²)/(2xy))²=-Delta²/(x²y²)`。
残り64fuv F²はdegree≤7、よってnormはdegree16、leadingは
`Delta^12/(4096*x^12*y^12)`で非零である。
大きい全base係数の新規展開は必要ない。

qR=qT=W=0なら、このnormからmはK0上代数的、さらにr,tもmonic quadraticから代数的。
するとxr+yt=k1/2-A0(m)がk1の超越性に反する。有限部にその共通零点は無い。

m=0の正規化上では、二節点の両枝でqRは単元（前の独立節点review）。
u=0またはv=0の場合はmがK0上代数的で、k2の独立性からm,f,他方のu/vは非零。
qR=0ならr=0、元Psiの式から `Psi=-m² f t²`。
t=0となりk1の超越性に反する。このidentityも独立に照合した。
f=0ではmuv≠0でqR=0からr=0、`Psi=uv(B0+2dt)²`。
d≠0かつtがk1について非定数なので零にならない。
c=0を分母として反転する必要はない。

残る有限零点ではmuvとqTが単元で、正規化のDVRにおいて
`qR=(W/m)²/qT`。従って零点の位数は全て偶数。
qRはaffine polynomialなので有限poleは無い。
K(m)上のfinite algebraで単元だったという旧結果を、projective曲線全体での単元性へ
取り違えていない。

## Infinity四点とRiemann–Hurwitz

旧probeのinfinity-positive部分は元のsmoothness-falseと区別して利用する。
保存された同じbinary quarticについて、独立な有限体Euclidean gcd(H,H')=1を確認する。
generic divisorはr⁴係数単元のfinite flat degree4 familyなので、generic infinityは相異なる四点。
`F(r,0,0)=16x⁴y⁸r⁴` は非零で、全四点のm homogeneous座標は非零。
intersection divisorがreducedなので、正規化上の各点でmは単純pole。

qRの全次数3部分は `Delta²/(4x²y²)*m³`。
r²のpole orderは高々2なので相殺できず、各点でord(qR)=-3。
節点解消はfinite部だけで、infinityのこれら四点を変えない。
従って二次cover h²=qRは四点でのみ分岐し、分岐次数は正確に4。
奇数poleを持つqRは幾何的関数体でも平方でなく、coverは幾何的にconnected/既約。

商genus1から `2g-2=2(2*1-2)+4`、ゆえにfull ordered curveの正規化genusは3。
既受理の全generic algebraはquartic×monic quadraticでdimension8。
今は両段が幾何的既約なのでこれはfieldとなる。
実際の選択像成分は同じ五baseへdominantで、内在的mも非定数であることが既受理の
universal六不変量独立性から従う。そのfield factorはこの唯一のfield全体なので、m-degreeは8。

結論はこのordered-base/relative-sign chartに限定する。
さらに係数体を代数的に拡大しても幾何的genusは変わらないが、Weyl groupで曲線を
新たに商したgenusを同一視しない。NEW multisectionの既知10/20/40、NEW定義式や
Poisson全成分尽くしとは別の結果である。

再現:

```
python code/q_firstjet_residual_genus3_branch_review_20261104_v63.py \
  --output /tmp/universal-ordered-genus3-review
```
