# 有限Sigma全域のraw四次rank78：split kernelと有限Nakayama

2026-11-03。既存のproper finite Sigma modelと、一回の実GF rank計算を組み合わせる。元の四quadの順序とa-clearing shiftsは `[0,2,2,1]`。Q2primeへのrow変更を混ぜない。

## 実測とsource binding

独立済みfinite-Sigma modelのp=32003、Sbarのdimension28、P/eのregular multiplication operatorsを再利用した。元24×56（四quad×六light）から四quadを剥がし、全24行を再構成して一致を確認した後、実84×126のquadratic-cofactor行列Mを作った。degree2/4 monomialはともにlexicographic exponent順である。

Mのregular representationは2352×3528、8,297,856成分。一回のFLINT rank計算は **2184=78×28** だった。rank時間2.279秒、全probe3.630秒。二回目rank・大きなinverse・別fieldは実行していない。全blockのu32 little-endian SHAはprobe receiptに保存した。compact inputsとnative codeから再生成できる。

## 六KoszulのsplitをQQで検証する

四quadのdegree2係数行列をQ（4×21）とする。既存の保存L（21×4）について、今回も全16のQQ恒等式

`Q*L=g*I4`, `g=−a²P²J`

を独自Fraction係数演算で再確認する。全2461cofactor項の数値分母はp-unitであることも検査する。Qの元804項を再構成し、modp後の全24×56がnative probe入力と一致することを確認する。

Kは六Koszulを行とする6×84行列である。行(i,j), i<jは、cofactor block iにq_j、block jに−q_iの21係数を置く。Mの行もblock i/cofactor monomialで同じ順序なので、K*M=0。これを全6×126のQQ係数恒等式として検査する。

Kのright inverse Rは、列(a,b), a<bについて **block aにL[:,b]/gだけを置く**。K*Rの(i,j),(a,b)成分はa=iならdelta(j,b)、a=jなら−delta(i,b)、他は零。順序i<j,a<bにより全体はI6となる。78角minorやArtin socleの推測を使わない。

## Proper finite V上の昇格

V=Z_(32003)、S=V[P,e]/(N,M0)を既存の有限modelとし、a=−P(20P+3e+20)/180を代入する。既存証明はe⁴/P²⁸のmonic relationsからSがfinite V-algebraであること、またFULLU=aPFB0CJがSのunitであることを示す。flatnessやgeneric fiber length28は仮定しない。

Q,Lの分母はp-unit、gはFULLUのunit因子の積なので、上のKとsplit RはS上で定義される。従って

`S^84 = ker(K) ⊕ S^6`

で、ker(K)はS上finite projective rank78、特にV上finiteである。このsplitによりmodp base changeはker(K)と可換となる。SのV-flatnessは不要である。

modpでdim ker(Kbar)=(84−6)×28=2184。観測されたMbarのrankも2184であり、Kbar*Mbar=0なので、image(Mbar)=ker(Kbar)である。よって

`C=coker(S^126→ker(K))`

はfinite V-moduleで、C/pC=0。有限NakayamaによりC=0となる。

QQへbase changeし、任意のcomplex Sigma点でKはrank6、Mはker(K)へ全射なので、**元raw84×126のrankは78** である。有限fieldでのrank下界を一般のQQ familyへ直接昇格したものではなく、同じproper finite modelとsplit kernelを介した有限Nakayamaである。

既存のfinite modelはFULLU-unitを全S上で与えるので、元ell chartとZW19もSigma上で適用できる。受理済みのdegree3 colon Koszul-lift argumentによりoriginal light rank146、weighted rank152もSigma全域で従う。

これは有限Sigma（G1の有限exception cover）に限る結果である。他の有限quad78例外、全aPJでのrank78、main C122のrank64、full quartic全射、長65、ambient除去は未証明のままである。旧omit-B²の撤回済みsyzygy推論は使わない。

再現（新しい出力先）:

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_finite_sigma_quartic_nakayama_20261103_v63.py \
  --output work/marked-square-v63/finite-Sigma-quartic-Nakayama-new
```

このdriverは保存GF rank receiptをbindし、QQ provenance/splitを検証するだけで、第二rankを起動しない。新しいGB・field・prime・大きなinverseは生成しない。
