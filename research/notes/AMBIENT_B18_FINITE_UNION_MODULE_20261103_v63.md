# generic B18のactual交叉：四行finite moduleとbounded target actions

2026年11月3日。前稿の共通R28・rank2 normalizer還元を用い、
実際のD=CA∩CBを二つの大きなimage idealのGBへ戻さず、四行の有限subquotientへ正確に還元する。
係数algebraic formal seriesの全展開やprimary構造は未計算である。

## 1. 実際の二成分とtarget

前稿により同じR=C[[r1,...,r28]]上に二つのfinite free rank2 algebraがある。

    BA=R[pA]/(pA²−tauA*pA−nuA),
    BB=R[pB]/(pB²−tauB*pB−nuB).

共通targetはT=R[[z1,...,z97]]で、各ziの二つの実際の像は

    zi^A=aAi+bAi*pA,  zi^B=aBi+bBi*pB.

全係数は前稿の有限29源circuit、28 unit implicit pivots、一つのquadratic Weierstrassから一意に定まる。
ここでsource normalizersをactual成分と同一視しない。
実際のimage algebrasをAA=image(T→BA)、AB=image(T→BB)とする。

B=BA⊕BBはR上free rank4である。各targetのpaired像をbi=(zi^A,zi^B)、diagonal unitを1=(1,1)と書く。
U=image(T→BA⊕BB)はactual unionの像、M=AA⊕ABとする。

## 2. polarized Cayley--Hamiltonによる有限word bound

RはQQ-algebraである。
rank4 free algebra Bでx=Σt_i biのCayley--Hamilton恒等式を取り、t1t2t3t4の係数を比較する。
先頭項は24*b1*b2*b3*b4であり、他の項はR係数と長さ3以下のproductだけになる。
24はunitなので、任意の四因子productは短いwordへ戻る。
長いwordの最初の四因子に反復すれば、全wordは長さ3以下で生成される。

係数を曖昧にしないため、S⊂{1,2,3,4}について

    c(S)=sum_partitions(pi of S) (-1)^(|S|−|pi|)
                       * product_blocks (|block|−1)! * Tr_B(product_block bi)

とすると、正確に

    24*b1*b2*b3*b4
      = sum_nonemptyS (-1)^(|S|+1)*(4−|S|)!*c(S)*product_remaining bi.

Traceはfree normalizer BのtraceなのでRに入る。
Uをfreeやflatと仮定していない。
この恒等式を二つのgeneric quadratic algebra上の四つのgeneric要素について、全四成分でsymbolically検証した。

従って97個のpaired target generatorsについて、Uは長さ0,...,3のwordだけでR-moduleとして生成される。
個数は

    1 + 97 + 4753 + 156849 = 161700.

全word indexを圧縮archiveへ保存した。各columnの係数は二次式での有限乗算だけで構成できる。
未知のstabilization stageを仮定しない。

同じrank2 polarizationは各Aiで

    2*b1*b2=Tr(b1)*b2+Tr(b2)*b1
                      −(Tr(b1)Tr(b2)−Tr(b1*b2))*1

を与えるので、AA,ABはそれぞれ1と97 target像の98個で生成される。
従ってMは196 columns、Uは161700 columnsの、同じ**四行R-matrix**で指定される。
rank2の二因子・三因子還元もuniversalに検証した。
Rはcomplete NoetherianでBはfinite freeなので、これら有限submodulesは閉である。
したがってformal power series像も同じ有限moduleとなり、polynomial像だけの主張ではない。

## 3. actual Dのexact module sequence

IA=ker(T→BA)、IB=ker(T→BB)は実image idealsである。
通常の二ideal exact sequenceから

    0 → T/(IA intersect IB) → T/IA ⊕ T/IB → T/(IA+IB) → 0

がある。最初の三つはU、AA⊕AB=M、Dそのものである。
よって

    underlying R-module of actual D = M/U

が正確である。最後の射は、各componentのtarget liftの差をIA+IBを法として取るものである。
非被約性やembedded部分を消していない。

重要な区別として、Uは一般にBのidealではない。
従ってB/Uをnormalizer componentwise multiplicationで環にする議論は使わない。
またnormalizer tensor productをactual Dと呼んでいない。

## 4. 有限target actionsが全環も復元する

ziのpaired multiplicationはBの基底(1,pA,1,pB)上で明示的な四次行列

    blockdiag( [[aAi,nuA*bAi],[bAi,aAi+tauA*bAi]],
               [[aBi,nuB*bBi],[bBi,aBi+tauB*bBi]] )

である。これらはMを保ち、Uも保つ。後者はword長4の上の恒等式で有限に検査できる。
よってM/Uに共通target Tの可換な作用が降りる。
unitに対応するclassは(1,0)∈AA⊕ABのclassである。
このclassはT-cyclicで、annihilatorは正確にIA+IBとなる。
したがって四行finite subquotientにこの97個のtarget actionとunitを付ければ、
actual完成環Dを全scheme structureのまま復元できる。
任意のnormalizer要素のmultiplicationを quotientへ降ろす必要はない。

全係数matrixの巨大展開を済ませたという主張ではない。
今回凍結したのは、前稿の明示有限source coefficient recipeから、
四行・有限個のword・有限個のtarget actionsだけでactual Dを復元する**有限かつ厳密な還元**である。
次は実際に係数を小さいnormalizer basisへ戻し、M/Uのnilradical・Ass・depthを計算する。

## 証明書

* `code/ambient_B18_finite_union_module_20261103_v63.py`
* `reports/AMBIENT_B18_FINITE_UNION_MODULE_20261103_v63/{report.json,finite-module-recipe.json,paired-words-degree0-through3.json.gz}`

universal symbolic検査1.39秒。
四因子rank4式、二因子・三因子rank2式、全161700 word indicesを保存した。
前稿の全source flags、共通parameter minor、quadratic coefficientの証明と組み合わせてgeneric B18に適用する。
full completed localtypeやprimary decompositionを解決済みとは呼ばない。
