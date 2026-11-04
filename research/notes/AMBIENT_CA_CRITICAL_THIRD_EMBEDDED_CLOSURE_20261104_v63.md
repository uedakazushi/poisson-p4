# CA-critical V 軌道は actual D の第三の相異なる embedded associated closure を与える

2026-11-04。D は被約成分 CA と CB の **実 scheme 交叉** である。
既存 V1 点の稠密 CA-critical 軌道を O_C、その閉包を Z_C と書く。
本稿は whole proper fiber と actual full contact completion の受理済み証明を
接続し、Z_C が既知 Pa、Qb の二つと異なる embedded associated closure である
ことを証明する。Poisson Jacobi scheme の新しい irreducible component や全
embedded-prime list の尽くしではなく、全成分下界58を変更しない。

## 1. 実在する既約20次元軌道と proper map の共通有限開

作用群を GL5 と非零 Poisson scalar の積とする。これは連結既約代数群である。
保存された同じ有理点 V1 の全 GL25 列は rank20 であり、scalar 列もその中に
入る。従って O_C は smooth irreducible な cone 次元20の locally closed 軌道。
その既約閉包 Z_C は20次元で、projective 版は19次元となる。これは新しい
点の rank 標本ではなく、既に独立監査済みの full150 source 証明を使う。

全 critical source の monic V 曲線は P0 から V1 を通り、その G 飽和閉包は
20次元、O_C はその稠密開軌道である。特に P0 belongs to Z_C。
この同定・両端・全 critical 方程式の零は CA V1 独立レビューで確認済み。

P0 における **全 proper** CA/CB source fiber は一台長2である。端点・全旗を
含む先行 Darboux 証明と全 fiber ideal を使い、選択 chart 内だけの数ではない。
properness により P0 を含む開で両写像は finite となる。その有限 module の
fiber 長さが2以下の locus も開であり、作用同変性によりこれらの開は G-stable。
Z_C は P0 を含み O_C が稠密なので、共通有限/長さ<=2開は O_C に会う。
G-stability により **軌道 O_C 全体がこの開に入る**。

この開上では両 proper source は実 normalization である。正則 source の
finite 性だけから birationality を仮定しない。P0 の受理済み completed image
A_i=R+I_i*tau_i では I_i!=0 なので tau_i=(f*tau_i)/f が image の分数体に
入り、finite source はそこで birational である。source と target の既約性から
generic degree1となり、有限開上の normal source が normalization を与える。
これは global proper map が至る所 finite normalization だという主張ではない。

## 2. CA は長さ2の一点、CB は滑らかな一点――両方とも whole fiber

CA source の V1 微分は corank1。したがってその選択点での局所 fiber は長さ
2以上であり、前節の全 fiber 長さ<=2から、全 CA fiber はその一点だけの長さ2。
全完成 Ac8 像の証明により fiber は正確に k[epsilon]/epsilon^2 である。
第二の遠い旗や compactified endpoint を残す余地はない。

CB については chosen smooth lift だけから whole fiber1と結論しない。
受理済みの common physical R0 conductor identity

    I_A0+I_B0 is maximal-ideal-primary

を使う。CA-critical の零weight V branch は、monic proper critical source と
common G19 cut によって Z_C の同じ branch と同定されている。その一般点では
I_B0 は単元。intrinsic conductor quotient は actual R0/I_B0 だから、これは
単なる抽象座標の一致ではなく、その一般 CA-critical 点が CB conductor の外
にあるという意味である。この source/branch binding は COMPARISON_C_GENERIC_DVR_RANK
の独立レビューで確認済みである。

より厳密に、もし Z_C 全体が共通有限開の CB conductor に入れば、その coherent
conductor ideal の完成と上記 actual V branch への制限も消えなければならない。
I_B0 が一般点で単元であることに矛盾する。従って O_C 上で conductor 外の開は
非空である。conductor と O_C は G-stable なので O_C 全体が conductor 外となる。
有限 normalization は conductor 外で target と同型であり、source は正則。
ゆえに CB は O_C 全体で smooth29、**全 CB fiber は一つの被約滑らかな点**。
既存の正確な CB29 differential はこの結論と一致するが、その rank だけを
whole-fiber 証明には使っていない。

全 generic pattern は次である。

| 稠密20次元 locus | 全 CA proper fiber | 全 CB proper fiber |
| --- | --- | --- |
| O_C | 一点、非被約長さ2 | 一点、被約長さ1、smooth |
| Pa の受理済み unit open と G 飽和 | 二つの被約 smooth 点 | 一つの被約 smooth 点 |
| Qb の受理済み unit open と G 飽和 | 一つの被約 smooth 点 | 二つの被約 smooth 点 |

Pa/Qb の全非零旗数 theorem と actual full-ring の非空 unit open を使い分ける。
full-ring 型を全非零 parameter へ拡張する主張はしていない。

## 3. full actual completion の associated prime は軌道閉包へ降りる

新しい source および純代数独立監査により、p in O_C で

    completed O_(D,p)
      = k[[g20,l7,v,x]]/(x*l7,x^2,x*v^2)

である。追加 ambient graph を全て残した全 ideal の等号であって、接錐だけ
ではない。非零 class x*v の annihilator は

    P_C=(l1,...,l7,v,x)。

従って P_C は completion の embedded associated prime。主 prime (x) の
商は smooth N28であり、P_C の商は regular20なので、これは minimal prime
ではない。monomial primary intersection は

    x*(l7,x,v^2) = (x) intersect (l7,x^2,v^2)。

作用群の20次元 section と ambient transverse slice の formal IFT product は、
CA、CB、N、その両 normalizer 比較を同時に保つ。その slice の原点を動かす
群 factor はまさに O_C の形式近傍である。従って P_C は代数的 locally closed
軌道 O_C の ideal の completion である。この同定により、一般の formal
fiber prime の contraction を無根拠に同じ次元と扱う問題を避ける。

faithfully flat completion の associated-prime contraction により
I_(Z_C,p) belongs to Ass O_(D,p)。localization の associated-prime 定理から、
global cone D の associated prime は I_(Z_C) そのものとなる。Z_C は既約20次元、
非零 homogeneous cone なので projective D にも embedded associated closure
Z_C^proj（次元19）を与える。

N は28次元であり Z_C はその proper subvariety。既知のもう一つの minimal
component R20とも同一視しない。実 N incidence により O_C subset N であり、
R intersect N=B18 の cone 次元は19。G-stable R が20次元軌道 O_C を含むことは
できない。この点からも Z_C が minimal R ではないことが分かる。

## 4. Pa、Qb の閉包との相違

既存の source morphism G times Pa_unit と G times Qb_unit の像は constructible
であり、その閉包 Z_P、Z_Q は既約20次元である。両 full D theorem は上表の
whole flag pattern を持つ非空 unit open とその飽和上で、embedded associated
prime を証明している。

仮に Z_C=Z_P とする。O_C と Pa_unit の G 飽和は、この同じ既約空間における
二つの dense constructible subset なので、それぞれ非空開を含む。二つの開は
交わるが、その点で全 CA fiber が一台長2かつ二つの被約点であることになり矛盾。
Z_C=Z_Q でも全 CB fiber の一点/二点が矛盾する。Z_P!=Z_Q は同じ既受理の
whole-fiber argument である。

これは fiber の台数の special limit だけを比較した議論ではない。一点長2が
二被約点の特殊化になること自体は認めるが、同じ20次元既約閉包の二つの dense
constructible locus が異なる generic pattern を持つことはできない。

よって actual D は、少なくとも三つの相異なる global embedded associated
closure Z_P、Z_Q、Z_C（cone20/projective19）を持つ。全 list の尽くし、閉包の
全境界、B18 原点の全 associated primes、Poisson 成分の追加数は主張しない。
新しい source/jet/点/GB/field 計算は行わず、exact依存とSHA256を同名 report に固定する。
