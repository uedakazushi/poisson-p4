# 平面二H＋cubicの残三全源：generic磁気でzero-surfaceが無くT38全親を除く

2026-11-08、PRIMARY E。純手のwhole-source必要条件。新operator、point評価、rank、minor、
profile、jet又はCASは無い。named陽性親又は全source完了はまだ主張しない。
旧source／T21陰性／whole126還元／global台帳／Gitは不変。

## 1. 全三lambda源とclosed zero-surface条件

全磁気六係数を失わないactual shear／逆は
`TWO_LINES_PLANE_CUBIC_THREE_MAGNETIC_QUOTIENTS.md` とその独立査読に固定する。
node-flex及びcusp-flexの全generic sourceはR=lambda h²+hy、
cusp接線・割線は3(lambda h²+y²)へ戻る。
元A6のJ又はI零補集合は別の同じsource閉境界として保持する。

projective Poisson zero schemeをZpiとする。
P4を使う普遍zero incidenceはproperなので、dimZpi≥2はparameter上のclosed条件。
真T38 sourceのgeneric全点にはraw zeroのtwisted-cubic cone surfaceがある。
実際 `reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json` でp0,t2,t3非零なら
\[
 x_1x_3=-2t_2x_4^2/p_0,\quad x_2x_4=2t_3x_3^2/p_0,
 \quad x_1x_2=-4t_2t_3x_3x_4/p_0^2
\]
上で全十entryが零。このP3内twisted cubicと自由x0のprojective coneはsurfaceである。
従って **T38全閉成分はdimZpi≥2** に入る。
以下で各全lambda線の稠密開がこのclosed条件の外にあると示す。

## 2. h非零：全anchorの共通零点を解く

h一のchartではh-Casimir liftがprojective bivectorを忠実に表す。
元全源P=U∧partial_z+V∧partial_w+(-d w+R)partial_z∧partial_wでは、
zeroはU=V=0かつ−dw+R=0と必要十分。
fullRの有限shearは全actualmapなので、この計算を元sourceにも戻せる。

**node-flex。** U,V,dはTHREE_MAGNETIC_QUOTIENTS§2の全式。
全base共通零点は [h:x:y]=[1:0:0],[1:1:0]のみ。
前者d=−1なので、fiberにはw=−lambdaという一本のprojective line。
後者d=0,R=lambdaなので、lambda=0の時だけ全fiber P2がzero。
lambda非零ならそのplane上のzeroは無い。

**cusp-flex。** 同note§3の全式から共通零点は
[1:-1:0] と [1:0:1],[1:0:-1] の三点のみ。
前者d=−1で一本のfiber line。
後二点はd=0でR=lambda±1だから、lambda=-1又は+1の対応する時だけ全fiber P2がzero。
それ以外のlambdaではその二planeにzeroは無い。

**cusp接線・割線。** 全式は
`TWO_LINES_CUSP_TANGENT_SECANT_MAGNETIC_NORMAL_FORM.md`。
h非零のbase共通零点は [1:0:1]だけ。
Pzw=3xz+3(lambda h²+y²)のx=0制限は3(lambda+1)。
従ってlambda=-1だけ全fiber P2がzeroで、他lambdaではh非零のzeroは無い。
別base点[0:0:1]は次節のh零側で扱い、除外しない。

## 3. h零：raw zeroとprojective zeroを区別する

h零かつx非零では、各全raw4×4行列のPfaffianはそれぞれ
−4x四乗、−6x四乗、9x四乗（いずれも非零定数倍で十分）。
従ってraw rank四であり、Euler wedgeで消えるprojective zeroにはならない。
rank二のdecomposable E∧Lへ等しくなることが不可能だからである。
単にh=0を捨てる議論ではない。

h=x=0,y非零では、node-flexとcusp-flexのraw tensorはそれぞれ
−3y²partial_y∧partial_w、+3y²partial_y∧partial_w。
Euler方向はy partial_y+z partial_z+w partial_wなので、
projective zeroであることとz=0は同値。この全zeroは一本のprojective lineでありsurfaceでない。

第三源ではraw tensorは3y²partial_z∧partial_w。
Eulerはy成分非零なのでそのimage planeに入らず、projectivezeroでない。
h=x=y=0のvertex P1は三源ともrawzeroで、全て保持する。
これ以外のh零点は無い。

## 4. Whole親への厳密な帰結

従って各稠密開
\[
 \lambda\ne0\quad(node),\qquad\lambda\notin\{-1,1\}\quad(cusp-flex),
 \qquad\lambda\ne-1\quad(cusp接線・割線)
\]
の全点でdimZpi=1。
元六磁気源のPGL飽和閉包SiがT38に全体として入れば、その全点がclosed条件dimZpi≥2を満たすはずで矛盾。
よって **三つのwhole23 sourceの親としてT38は不可能**。
上記稠密開の個別lambda点もT38外である。

例外lambda=0／±1／-1にはzero planeが実在し、T38所属も非所属も本稿から主張しない。
元J／I零境界にもgeneric陰性を外挿しない。
全sourceを既知有限閉和へ入れるという陽性親認識は依然別義務である。
既 `TWO_LINES_PLANE_CUBIC_KNOWN_PARENT_REDUCTION.md` のwhole126還元を使う場合、
残known24親はT36又はT37へ縮むが、whole現component集合が増えた後の比較を自動更新したとはしない。
新しい成分性や全E完了は結論しない。
