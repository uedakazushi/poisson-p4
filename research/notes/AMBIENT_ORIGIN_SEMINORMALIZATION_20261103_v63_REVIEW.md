# 全形式CA・CBモデルのseminormalization：独立検証

2026-11-03。受理済みの全形式モデルを前提に、`AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63`
を承認する。個々の原点completed local ringに限り、global/relative embeddingの分類ではない。
producerは再実行せず、独立の整数even/odd pair算術と一般のsquare/cube criterionを用いる。

R、B=R[t]/(t²−nu)、Iは全モデルの通りとし、
`F=nu−s(s−r)`、`J=(z6,F,w,r)` とする。
Jはself supportのprimeで `R/J=C[[g19,s]]`、nu=s²。
radical(I)全体とは異なる。I+(w,r)=Jなので、
`R directsum J*t=A[wt,rt]` はalgebraとして正確に生成される。

独立word算術で二つのcube恒等式を全係数確認する。
`nu*w³=F*w³+(ws)w²(s−r)`、`nu*r³=F*r³+(sr²)r(s−r)` により、
wt、rtのsquare/cubeは両方とも元Aに入る。
b²,b³を持つelementのadjoinはelementary subintegralである。
実際各residue fieldで、b²≠0ならb=b³/b²、b²=0なら唯一のresidue値b=0となり、
integral mapはprime上全単射かつresidue field不変。
従ってAから候補への二段階拡大はfinite subintegral。
w/rは各critical conductor軸で非零なのでIに入らず、拡大は真である。

## normalization全体でのseminormality

任意のb=x+yt∈Bに対し、square/cubeのodd係数は正確に
`2xy`、`y(3x²+nu*y²)`。b²,b³が候補に入るなら、これらはJに入る。
domain `R/J=C[[g19,s]]` へ落とし、y_bar≠0と仮定する。
第一式はx_bar=0を強制し、第二式は `s²*y_bar³=0` となって矛盾する。
sはunitでなくてもdomain内の非零元であり、この議論は原点でも有効。
従ってy∈J、bは候補に入る。

total fraction fieldの任意の候補でsquareが候補ringに入れば、そのelementはintegralなので、
共通normalization Bに入る。従って上の計算はB内の一部だけでなく全square/cube criterionを覆う。
候補ringはreduced domainでchar0、ゆえにseminormal。
finite subintegralなseminormal intermediate ringなので、正確にAのseminormalizationである。
1∉Jなのでtは入らず、候補はBより真に小さい。

## conductor・全81 kernel・不変量

候補のconductorはtargetで `J directsum J*t`、sourceでJB。
BではF=t²−s²+srでありr∈Jなので、
`JB=(z6,w,r,t²−s²)`。
source conductor quotientは `C[[g19,s,t]]/(t²−s²)` という被約node。
target quotientはsmooth `C[[g19,s]]` で、generic二点gluingは保持される。
これはBそのものをsubintegral normalizationとするNの別定理とは異なる。

Fへのcoordinate change後、z6,F,w,rは正則九生成列。
独立コードはJacobian determinant1を確認する。
全第一syzygyは36 Koszul、pair productsは45。
integer pair算術で保存された全81関係を照合し、全てが零になることを確認する。
正則列のKoszul完全性と、55の場合と同じ収束するR+sum RY還元により全kernelとなる。
GBやM2は不要。

depth J=21、従って候補ringのdepth21。mu(J)=9なのでedim29+9=38。
odd productがm_R³に入り `m_candidate²=m_R*m_candidate`、従ってm_Rはreductionとなる。
R-rank2を用いてe=2。cone dim29/edim38/depth21/e2、projective dim28/edim37/depth20/e2。
元のCA/CB座標への相対alignmentやglobal seminormality locusは主張しない。

独立コードと証拠：`code/ambient_origin_seminormalization_review_20261103_v63.py`、
`reports/AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63_REVIEW/`。
