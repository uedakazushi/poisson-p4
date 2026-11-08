# 残六候補の相互包含：proper divisor incidenceと既約次元による有限監査

2026-11-08、PRIMARY exceptional membership／component-parent audit。
**producer純監査、独立review待ち。** 新point／operator／minor／jet／solver無し。
rootだけがcatalogue/Git/mainを統合する。本稿は候補をcomponentとして加算しない。

## 1. 対象と既存下界・上界

Kを凍結prior142の閉成分の和とする。旧
`KNOWN142_SEVEN_CANDIDATE_WHOLE_POISSON_COVER_AUDIT.md`はP_red=K∪七既約候補を証明した。
その旧149上界proofは変更しない。今回wholeh5の残C_UはT19へ入り、Kへ吸収される。
従って残る候補は以下の六つだけである。

|記号|正確source|dim|
|---|---|---:|
|C_G|Ga9 linear-w unit canonical actual orbitの閉包|23|
|C_N|node_flex_tangent_through_flexの全A6磁気源のPGL飽和閉包|23|
|C_F|cusp_flex_secant_through_flexの全A6磁気源のPGL飽和閉包|23|
|C_C|cusp_cusp_tangent_cusp_secantの全A6磁気源のPGL飽和閉包|23|
|C_D|constant DOUBLEの全sourceを正規化した一actual orbitの閉包|23|
|C_T|constant TRIPLEの全sourceを正規化した一actual orbitの閉包|22|

すべて既約である。actual orbit閉包は連結PGLの像閉包、plane三つはPGL×A6の全域像閉包。
Ga9は既独立M48/Euler24でorbit23、全11源/磁気unit正常化により全像閉包がこれと一致する。
plane三源は各fullsourceの実方向rank48が与える下界23と、
固定D fiber≤9＋base/essential-frame≤14の全上界23が一致する。
constant DOUBLEは既全51行列M48からorbit23、TRIPLEはM47からorbit22。
これらは既保存同一点行列のrank下界の再使用で、新点を一般点と見做した下界ではない。
候補が全Poisson成分であるという主張はここからは出ない。

## 2. generic actualDの全source不変量

plane三源の座標は(X,Y,Z)、K3=Y²Z−X³、N3=Y²Z−X²(X+Z)とする。
全R6を動かしてもsourceのDは変わらない。

|候補|dense sourceのactual divisor（scalar除く）|reduced|本質変数数|
|---|---|---|---:|
|C_G|h(4h³w−8h²y²+y⁴)|yes|3|
|C_N|Z(X+Z)N3|yes|3|
|C_F|Z(X−Z)K3|yes|3|
|C_C|Y(Y−X)K3|yes|3|
|C_D|h²(x²y+h²x+h³)|no|3|
|C_T|h²(x³+3hxy−3h²y/2+h²r)|no|4|

Ga9 quarticはprimitiveなw一次式、DOUBLE cubicはprimitiveなy一次式、
TRIPLE cubicはprimitiveなr一次式なのでそれぞれ絶対既約。
plane N3/K3はそれぞれirreducible nodal/cuspidal cubic。
各五次式は表示したessential subspaceを本当に使う。
TRIPLEについてはDのconstant directional derivativeを零とすると、
r項からh方向、y項からx方向、x項からy方向、最後にr方向が順に零となる。
したがってDのessential rankは正確に四であり、三以下ではない。

C_FとC_Cは両方cuspidal cubic＋二lineであるが、unordered factor pairの配置が違う。
K3の唯一cuspはO=[0:0:1]。
C_FのZとX−ZはどちらもOを通らず、交点I=[0:1:0]はsmooth flex。
C_CのYとY−Xは両方Oを通り、前者はcusp tangentである。
従ってこの二つのactualD配置はPGLで同じではない。

## 3. 全boundaryを含むproper closed条件

V*=H0(P4,O1)とする。以下の像はprojective parameter spaceからのproper像なのでclosedである。

* M_113={L1 L2 F3}：二linear factors＋任意cubic。L1=L2、F3の分裂・zero-residueを除かない。
* M_23={L² F3}：double hyperplane＋任意cubic。F3自身にL因子がある境界も含める。
* Sigma3={D∈Sym5(W),dimW≤3}：Gr(3,V*)上のprojective bundleをproperに忘れる像。

actual Pfaffian mapはrank四open（D非零）でregularなので、その逆像もこのopenでclosed。
rank二boundaryではD0であり、affine coneの零を加えて解釈すれば矛盾しない。
rank二境界の有無を本稿の非包含から排除しない。

C_N/C_F/C_C全体のrank四部分はM_113∩Sigma3、
C_D/C_T全体のrank四部分はM_23、C_DはさらにSigma3に入る。
C_Gも全rank四部分がSigma3に入る。
したがってC_GのdenseDはM_113の外、四つのreduced denseDはM_23の外。
一方double-H cubicはM_113にも入り得るし、1+1+3はM_14にも入る。
この片方向を逆にして因子分割だけで誤った排除をしない。

さらにcuspidal-cubic／二lineの正確flagを使うと、
「mark qでcubic singular、tangent cone rank≤1、両lineがqを通る」incidenceもproper closed。
C_Cはこの条件内だがC_Fは唯一cuspが二lineに乗らずその外。
node N3は唯一singularityがrank二tangent coneなので、cuspidal-cubic incidenceの外。
これらのflagはsimple boundaryだけを調べた非包含ではなく、印刷式で決まるwhole closed条件である。

## 4. 同dimensionの五候補は全て相互非包含

C_G,C_N,C_F,C_C,C_Dは既約dim23。
もしC_i⊂C_jならdim一致でC_i=C_j。
それぞれのrank四sourceのD-imageは上表の一つのPGL orbitにあり、
このsourceはC_iにdenseである。したがって等しい閉包なら、そのgenericD配置も一致しなければならない。

上表で異なるものは、prime-degree/multiplicity partition、node/cusp、
又はcuspidal cubicと二factor linesのcusp incidenceのいずれかで区別される。
従って五つは全て相異なり、どちら向きの全閉包包含も無い。
構成sourceの一つの特殊Dだけをgeneric化したのではなく、全磁気sourceで同一Dという証明を用いている。

## 5. 22次元TRIPLEと他五つ

C_TのdenseDはessential rank四。
他五候補のwhole rank四部分はproper closed Sigma3内だから、
C_Tのcanonical rank四点はどれにも入らない。従ってC_T⊄C_j（j≠T）。
逆にdim23>dim22なのでどのC_jもC_Tへ入らない。
ここはgenericDのessential-variable下界とwhole parentのclosed上界の正しい方向である。

## 6. 最小相互包含表と分類への意味

行が『含まれる候補』、列が『候補親』。○は同一、×は全閉包包含無し。

|⊂?|C_G|C_N|C_F|C_C|C_D|C_T|
|---|---|---|---|---|---|---|
|C_G|○|×|×|×|×|×|
|C_N|×|○|×|×|×|×|
|C_F|×|×|○|×|×|×|
|C_C|×|×|×|○|×|×|
|C_D|×|×|×|×|○|×|
|C_T|×|×|×|×|×|○|

この表は**候補同士が交わらない**という意味ではない。
magnetic-zero、factor衝突、rank二、他の低次boundaryは共通部分を持ち得る。
sourceごとの全点に親非所属を外挿せず、dense source／whole closure inclusionの判定としている。

142+六の有限被覆とこの表を合わせると、未知の成分数は

```
142 + #{i : C_i が既知閉和Kに含まれない}.
```

実際、C_iがKに含まれれば有限既知成分のいずれか一つに含まれる。
含まれなければ、有限被覆・既約性からC_iはある被覆項へ入り、相互非包含によりC_i自身が全成分となる。
従って六つのknown-parent membership判断は相互吸収で代用出来ない。
本稿はその六つの判断自体、T37/T39等への陽性、追加成分性、最終exact総数を結論しない。
旧149被覆audit・catalogue・mainを改変せず、この最小posetだけをrootへ渡す。
