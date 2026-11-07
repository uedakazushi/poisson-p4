# normalな次数六rank零曲面に残る幾何：最大sectional genusと非有理特異点の排除

2026-11-07、rootによる幾何証明。主義務Gの部分還元であり、全次数六曲面の分類ではない。
別作者の確認前にはglobal主張台帳へ昇格しない。

## 1. 命題

S⊂P4をnormal・integral・非退化な次数六曲面とする。
既約被約五次の実退化因子を持つrank四Poisson構造πがS全体で零なら、

1. q(S)=h1(O_S)>0。
2. 一般hyperplane曲線Cの種数gは三以下。
3. Sの特異点は全てrationalで、minimal resolutionのKodaira次元は−∞。
   その解消は正種数q(S)の曲線上のruled surfaceのblowupである。

横断node、clean交叉、正規化又は導手のsmoothnessを仮定しない。
この残るirregular ruled解消を持つnormal surface自体の不存在は、ここでは証明しない。

## 2. sectional genus四ならCI(2,3)でありq=0

normal surfaceの特異点は有限個なので、一般hyperplane曲線Cは滑らかで非退化P3曲線となる。
degree六のCastelnuovo上界はg≤4である。
g=4と仮定する。L=O_C(1)はdegree六、Riemann–Rochは

    h0(L)=3+h0(K_C−L)

を与える。h0(L)≥4、deg(K_C−L)=0からL≅K_C、h0(L)=4となる。
特にh0(O_C(2))=12−4+1=9。

INTEGRAL_SURFACE_LINEAR_NORMALITY.mdにより、degree六のSはlinearly normalでh0(O_S(1))=5。
完全列を用いて

    h0(O_S(2)) ≤ h0(O_S(1))+h0(O_C(2)) ≤14

となるから、Sは非零quadric Qに入る。
Sのintegral非退化性からQはirreducible、従ってrankは3,4,5のいずれか。

rank五ではPic(Q)=Z H、Sのdegree六から[S]=3H。
rank四ではCl(Q)=Z[P1]⊕Z[P2]、H=P1+P2。
一般hyperplane上のsmooth quadricへの制限は、Sのclass(a,b)をCのbidegree(a,b)へ送る。
degreeとgenusから

    a+b=6, (a−1)(b−1)=4

なのでa=b=3、やはり[S]=3H。
rank三ではCl(Q)=Z[P]、H=2P、degree六から[S]=6P=3Hである。
この二つのsingular quadric class groupとrulingの記述は、既に一次照合した
Tironi arXiv:1201.4508v1、printed p.5のquadric coneの計算と同じものである。

いずれもSはCartier divisor3HとしてQのcubic sectionである。
H0(P4,O(3))→H0(Q,O_Q(3))は全射なので、Sはscheme-theoretic CI(2,3)となる。
Koszul完全列からq(S)=0を得る。
ここではcycleの次数だけでlinkage又はnormalizationを同一視していない。

## 3. rank零quintic仮定はq=0を全て除く

degree六・q=0なら、Euler全source補題とlinear-normalityにより
H0(Λ2T I_S)=I_S(2)⊗Λ2 C5。
二つの独立quadricはSを含めないため、全πの実Dは零又はq2²linearとなる。
これは既約被約Dの仮定に反する。
従ってq(S)>0。§2からg=4も不可能なのでg≤3。

## 4. 非有理特異点はSchröerの定理に矛盾する

H=O_S(1)とし、normal canonical Weil classとHとの交点を取る。
CはSregに入り、通常のcurve adjunctionから

    K_S H=2g−2−6=2g−8≤−2。

Hがampleなので、全n>0でnK_Sはeffectiveにならない。
nonrational singularityがあれば、S. Schröer,
*Normal del Pezzo surfaces containing a nonrational singularity*,
Theorem2.2（printed pp.4–5）の一次確認済み定理によりq(S)=0となり矛盾する。
この定理はその名前とは異なり、当該statementにはdel Pezzo又はPicard rank一の仮定が無い。
proper normal surface、perfect field、nonrational singularity、全nK ineffectiveだけを用いる。

従ってSの全特異点はrationalである。
minimal resolutionf:S̃→Sについてq(S̃)=q(S)>0。
pullback f*Hはnefかつbigで、K_S̃ f*H=K_S H<0なので、全pluricanonical classはineffective。
よってκ(S̃)=−∞、標準的な滑らか曲面のruledness定理から主張3を得る。

## 5. 未確認文献と残る課題

E. Halanay, *Normal surfaces of degree 6 in Pn*,
Rev. Roumaine Math. Pures Appl.34 (1989),111–116の書誌は、
Halanayの別論文のpublisher referenceで確認できる。
本文はまだ取得・精読できておらず、その分類結論は上の証明に一切使っていない。
調査したpublic archiveは2007年以降しか示さず、1989年本文への入口は未発見。
non-normal degree5/6論文 <https://doi.org/10.1007/BF00147920> のpublic PDFも取得できなかった。
具体的にはweb readerがpublisher PDFを取得できず、既設Cloud proxyを使った通常curlでも
`https://link.springer.com/content/pdf/10.1007/BF00147920.pdf` へのCONNECTが403で拒否された。
TLS trustは変更せず、非公開copyや別経路によるアクセス制限回避は行っていない。
これらの本文を正規の公開copy又は認証済み購読経路で読めれば、残るruled型の
全normal-model分類を照合できる可能性がある。abstract/referenceを証明扱いしない。

未処理はrational singularitiesを持つirregular ruled解消のnormal degree六モデルである。
滑らかな場合は既存の小次数irregular曲面補題で次数五楕円scroll以外を排除済みだが、
そのsmooth double-point formulaを特異Sへ無断適用しない。
非normal S及びrank二restrictionも本noteの射程外に残す。
