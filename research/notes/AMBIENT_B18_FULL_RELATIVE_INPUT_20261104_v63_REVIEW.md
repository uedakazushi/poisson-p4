# Actual full186入力と四加群行列：独立レビュー

2026-11-04。**CBの完全な99 kernel生成元、全87 actual graph差、四Q[n]とJinfの
有限生成行列を受理する。** literal column countをrankや自由性と取り違えない。
所属・colon・socle・depth19/20・nilradicalの判定はこのpacketに含まれない。

## Actual座標と非零g

前stageのactual source/48unit系と独立TOP12レビューをhashで固定して再利用する。
F9=(zero7,z; j46,j95)は同じphysical R0上の正則座標でdet=1/96。
ここでは(z1,...,z7,r,s)と書く。除いたf63は一次部分が零。
pointed inverseをz=0へ制限したg(r,s)と、七coordinate telescopingから
`f63=g+Σz_i e_i` を得る。全16 inverse instanceの全912式と変数bindingをread-onlyで再構成し、
全784 coordinate quotientの分子・正則座標・label・原点値を照合する。
新しいsource、author builder、GBは再実行しない。

g=0ならactual ideal I_CBはz7,wr,wsの九本で生成される。
独立に受理済みのactual μ(I_CB)=10に反するためg≠0。
これはgを有限jetから非零と推測したものではない。
gはr,sの函数で、g=rA+sBの二つの正則coordinate quotientも明示される。
f63の一次部分零からA(0)=B(0)=0という保存値も正当である。

## 四つのgによる除算は循環していない

共通minimal20のphysical functionsは各actual完成像のcotangent生成系であり、
形式的Cohen/Nakayamaによって全完成像代数を生成する。これは既存のactual埋め込みの結果。
非負weightのためweight0部分を生成するのはR0と八つのzero odd functionsだけである。
積のodd係数はその八係数のideal I0に属し、形式極限でもNoetherian idealが閉である。
したがって四zero graph出力のodd係数bは、未知の相対liftを仮定せずI0に属する。

z=0へ制限するとbbar∈(g)。g≠0かつR0/(z7)=QQ[[r,s]]がdomainなので、
RegDiv(bbar,g)は一意な正則power series。
分子・gが代数的なため商も代数的である。
これで具体的な四商函数が指定されており、未指定hを立てたmembership問題とは異なる。
そのorigin値は未評価のままでよい。

checkerは四g nodeをNoneとして保持し、加減乗算にのみ伝播させる。
全inverse nodeについてg ancestorが無く、origin値が既知の非零単元であることを確認する。
未評価の値を零にしたり、gまたはg商を単元と仮定していない。
残り83 graphはF9 coordinate telescopingで持ち上がり、全87 liftとgraph pairを保存nodeへ一致させる。

## 全44 first syzygyと55 productの完全性

実generator changeはunipotentで、actual idealは
`I=(z1,...,z7,g(r,s),wr,ws)`。
coreの関係 `a*g+b*wr+c*ws=0` をmod wで見ると、g≠0で商環がdomainなのでw|a。
a=w*hとしてwを消し、g=rA+sBを代入すると
`r(b+hA)+s(c+hB)=0`。
(r,s)が正則列なので、このsyzygyは

```
(w,-A,-B), (0,s,-r)
```

の二本で尽くされる。HB matrixのunsigned minorsは(g,-wr,ws)、signed minorsは
実generator順(g,wr,ws)となる。符号を取り違えない。

z7を戻す際、任意syzygyをmod(z7)で二core関係に還元し、liftして引く。
残るcore係数のz成分を7×3のcross関係で移し、最後のz7だけのsyzygyは
正則列z7の21 Koszul関係で尽くされる。
したがって21+21+2=44は完全なfirst-syzygy生成系である。
抽象free kernelのrankだけから全点での完全性を推論していない。

奇変数10個の全55 pair productsを加えると、任意式はscalar+linear oddまで還元できる。
自由basis(1,tau_B)への像が零ならscalar係数は零、linear係数は上の44 syzygyで尽くされる。
よって99本はCB完成像のkernel全体を生成する。99の最小性は不要。

## Same physical targetへの代入と四行列

CB odd targetへ代入するのは、CAの同じphysical出力からCB even partを引いたpairである。
f63のunipotent変更もこのpairへ施す。
奇pairの積はtau_A²=nu_Aを使い、CB product式のscalar右辺はnu_Bを使う。
checkerは両者を区別して全99式を再構成する。

全87 graph差ではweight0四本も残し、positive53/24/6本は適切なw冪でhomogeneousにliftする。
全入力はCB/CA各trace coverを同一視せず同じphysical base上で定義される。
graph lift選択差はCB kernelに吸収され、全kernelとlower graphを保持するので
actual D=A_CA/Qを失わない。

kernel degree分布64/32/3とgraph分布4/53/24/6から、Qの全186生成元は68/85/27/6。
checkerは全pair output IDを照合する。
CA C0の九生成元、Hの十一生成元を掛けて、各R0² columnも全て再構成する。

```
Q[0]:9*68=612
Q[1]:11*68+9*85=1513
Q[2]:11*(68+85)+9*27=1926
Q[3]:11*(68+85+27)+9*6=2034
Jinf:11*186=2046
```

全8,131 columnsと追加top12列を保存DAGへbindする。これらはliteral生成column数であり、
module rank、独立性、free性、所属、socleは計算していない。
以前の四加群depth criterionに必要なactual有限入力が揃ったという結論までを受理する。

再現:

```
python code/ambient_B18_full_relative_input_review_20261104_v63.py \
  --output /tmp/B18-full186-independent-review
```
