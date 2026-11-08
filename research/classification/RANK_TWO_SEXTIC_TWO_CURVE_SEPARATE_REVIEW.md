# Degree六rank二支持の全空性：二曲線留数証明の独立読解

2026-11-08。PRIMARY E。**純証明PASS**。
対象はrootの `RANK_TWO_SEXTIC_TWO_CURVE_RESIDUE_EXCLUSION.md`
（SHA578e49dc…）。新CAS、射影、point、operator、minor、rank又はjetを実行していない。

被約既約actual quinticを持つrank四Poisson構造について、
degree六integral非正規支持の一般restriction rank二という全枝は空となる。
旧G上界をpointwise Eへ読み替えるのではなく、旧全源が満たす留数必要条件へ
独立な二本目の曲線条件を加える証明である。

## 1. 使う先行全源定理と、その補集合

`RANK_TWO_SEXTIC_SURFACE_SEPARATE_REVIEW.md` §§2–4のrank非依存入口を使う。
degree六はδ_S一、section genus零と一つだけの二次元非正規支持Sを強制する。
正規化の完全系モデルは113、122、014、023、005の全五scroll型である。
005はline vertexとplane curveのjoinで、二次元非正規支持はplaneの有限和だけなので除かれる。

他四型をsmooth Y=P_P1(O(e0)⊕O(e1)⊕O(e2))上で扱う。
014/023はzero-section曲線を縮めるsmall resolutionであり、divisorは縮めない。
Pic(Y)=ZL⊕ZF、Σei=5、C~−KY=3L−3F。
generic transverseδ一の全型はnode/cuspで、cuspならC=2ΓとなるがL係数三の奇数性に反する。
nodeのconductorはprime一つ又はprime二つだけである。
有限正規化はdivisorをcurveへ縮めないため、generic Sを覆わない余分なconductor divisorを隠していない。
connected non-split node coverもprime一つの枝へ入る。

上記モデル尽くし、generic nodeでのPoisson lift・tangency、閉inverse-residue形式αの存在は
先行独立全源定理への明示的依存であり、今回再計算したとは主張しない。

## 2. αに必要なのはgeneric log性だけ

各conductor prime Γのgeneric smooth locusで
`α=λ dz/z+regular`、λ非零constantである。
restriction rank二とΓ→Sのgeneric separabilityがλ非零を与え、closednessが定数性を与える。
係数体をPoisson中心とは置かない。

Yの他のcodimension一にはpoleがない。
small mapはexceptional divisorを持たないので、pullbackで追加のdivisorial poleは生じない。
全域log-smoothnessを仮定する必要はない。
generic log normal formの失敗点、Γのsingular locus、Γ同士の交差、その他の未制御点は
有限個のcodimension二以下の閉集合へまとめ、それらを曲線から避ければ十分である。

## 3. 二種類のcomplete curveは実際に全bad locusを避けられる

Lはbasepoint-freeである。一般二hyperplaneのpullbackの交わりBHはsmooth curveになり、
bad curveの像を避ける一般codimension二線形空間を選べる。
small vertexは一点なので、それも避けられる。
BHが非連結なら各成分のresidue theoremを足せば同じ関係を得る。

一方、一般ruling fiber Yt=P2に各bad curveが残すのは有限集合である。
baseを支配しないbad curveは有限個のbase値だけを占めるため、一般tで避ける。
small exceptional zero-sectionもYt上の一点だけである。
その有限集合を避け、各Γi|Ytのgeneric smooth locusを横断する一般line Bfを選べる。
char零のgeneric smoothnessとBertiniにより、非被約なgeneric intersectionを仮定せず選べる。

したがってαの両曲線へのpullbackは、指定Γiとの交点以外にpoleを持たない。
各横断交点でresidueはλiそのものなので、完全曲線の留数定理は
`Σ(B·Γi)λi=0`を与える。射影ruling lineの無限遠点も含む全P1上の関係である。

## 4. 一primeとsplit二primeを全て排除する

C=ΓならL²Γ=L²(3L−3F)=15−3=12。
BH上の留数和が12λ零を与え、非零性に反する。

C=Γ1+Γ2なら、各ΓiはSへdegree一で写るのでL²Γi=6。
BH上でλ1+λ2=0。
Γi=aiL+biFと書くと、aiは整数でai=deg(Γi|Yt)、a1+a2=3。
Bf上ではa1λ1+a2λ2=0である。
二関係の行列式は6(a2−a1)で、和が奇数三なので非零。
従って両λが零となり、generic rank二に反する。

各aiの個別分類、|2L−4F|のprime存在、射影の一般位置又は導手sectionの具体的因数分解は不要。
旧023の全split-conductor proper locusにも同じ二関係を適用できる。
smooth支持Sを追加に仮定せず、node／cusp／non-split／split／全finite projectionと
codimension二導手境界を先行全入口から戻している。

以上で指定degree六rank二の全枝は空であり、他degree、rank零、可約actual D又は全研究の完了を主張しない。
旧G証明書・witnessは不変であり、新しい証明はそれらの数学的誤りを主張しない。
pure jobのgateは入力hashと契約だけを認証し、この全数学読解と区別する。
