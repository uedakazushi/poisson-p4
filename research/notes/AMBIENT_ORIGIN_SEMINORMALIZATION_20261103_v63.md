# full formal CA・CB modelのseminormalization

2026-11-03。前提は新しい全形式モデル

    R=C[[g19,z6,nu,w,r,s]],  B=R[t]/(t²−nu),
    I=(z6,F=nu−s(s−r),wr,ws,sr²),  A=R directsum I*t

である。actual CA・CBへの適用は、その全形式モデルの独立監査を通して行う。
globalやrelative CA/CB embeddingについてのseminormality分類ではない。

## 1. 保持すべきself branch

J=(z6,F,w,r)と置く。これはR/J=C[[g19,s]]、nu=s²という**self supportのprime**である。
全conductor idealのradicalを代わりに用いると、critical枝で必要なsubintegral補正が残るので正しくない。

候補は

    A_sn=R directsum J*t=A[w*t,r*t]

である。I+(w,r)=Jなので、この二つの新しいelementで正確にA_snが生成される。

## 2. Aからはelementary subintegral

(wt)²、(rt)²はRに入る。さらに正確な恒等式

    nu*w³=F*w³+(ws)*w²*(s−r),
    nu*r³=F*r³+(sr²)*r*(s−r)

があるので、(wt)³、(rt)³もAに入る。
従ってA→A[wt]→A_snはelementary subintegralな有限拡大である。
一方w,rはC=R/IのそれぞれW-axis、V-axis上で非零なので、wt,rtは元のAには入らない。

## 3. B全体でのsquare/cube test

b=x+yt∈Bとする。b²,b³がA_snに入るなら

    2xy∈J,       y(3x²+nu*y²)∈J

である。domain R/J=C[[g19,s]]へ落とす。
もしy_bar!=0なら第一式からx_bar=0、第二式からs²*y_bar³=0となり矛盾する。
従ってy∈J、すなわちb∈A_snである。
sをunitだと仮定していない。domain上で非零因子であることだけを使う。

total fraction fieldの任意elementでsquareがA_snに入るものはintegralなのでnormalization Bに入る。
従って上の検査はBの一部だけではなく、seminormalityのsquare/cube criterion全体を満たす。
A_snはseminormalであり、Aからの有限subintegral拡大でもあるので、これが正確にAのseminormalizationである。

## 4. normalizationとは異なる

tはA_snに入らないので、A_snはBより真に小さい。
source conductorは

    JB=(z6,w,r,t²−s²)B

で、target conductorはJ directsum J*tである。
target conductor quotientはC[[g19,s]]、source conductor coverは

    C[[g19,s,t]]/(t²−s²)

という被約nodeである。
従ってgeneric self branchの二点gluingは残る。
Nの以前のsubintegral normalizationと異なり、ここでBそのものをseminormalizationとは呼ばない。

Jは正則九生成列なので、A_snの全kernelは36 Koszul関係と45 pair products、計81関係で与えられる。
coneではdim29/edim38/depth21/e2、projectiveではdim28/edim37/depth20/e2である。
全形式モデルの座標g19のうちscalar方向を除くため、projectiveではg18となる。
新しい関係をGBで推測したものではなく、regular-sequence syzygyとR-basis(1,t)で完全性を証明する。
