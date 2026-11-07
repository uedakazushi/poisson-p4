# tangent conic・h=L^3：一つの元Jacobi式で全essential三へ

2026-11-07。主義務G/E。実D既約被約quintic、rank二restriction plane、
reduced tangent-conic zero cubic、normal位数二、constant double rootのh=L^3型を扱う。
一回の有限jobで唯一粗上界27となった全sourceである。
既complete K1と全九P2の支持、及び元projective Jacobiの一式だけで認識する。
新profile/minor/GB/jet又はspectral計算は行わない。

**定理。** 指定全fullhigher sourceのactualDは三つのlinear formsだけに依存する。
irreducible rank四部分は既essential-three quinticの八23次元成分へ全pointwiseに被覆され、G像は23以下。
可約又は非被約actualDの全E被覆を追加しない。

## 1. 全sourceのOre形

座標(x,y,z,u,v)、f=z(xz-y^2)、N=2y partial_x+z partial_yとする。
P0はX=-8N、Y=-5N、Z=z partial_x、normal bracket零。
`reports/CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_20261107.json` のconstant_hcubeにある
全350x30 kernel九次元をn0,...,n8で戻す。
P2はinhomogeneous weight二式を解く必要があるため、ker J2三次元だけに狭めず**全九P2**を許す。

W=(z,u,v)、B=C[z,u,v]と置くと、この全必要sourceは正確に

\[
 \{z,u\}=\{z,v\}=\{u,v\}=0,
 \quad\{y,W\}=\delta(W),\quad\{x,W\}=y A(W)+V(W)\tag{1}
\]

という形を持つ。delta,VはB上のquadratic coefficient vector fields、Aはlinear。
全K1/P2の支持からx又はyへの他の依存は存在しない。
具体的に

\[
 A=\bigl(2z+(9n_6+3n_8)u+(2n_2+4n_4)v\bigr)\partial_z
       -16u\partial_u-10v\partial_v,
\]

\[
 \delta_z=z^2+n_2zv+\tfrac92(n_6+n_8)zu+P_{2,yz}(u,v),
 \quad\delta_u=-8zu+n_4uv,\quad\delta_v=-5zv+n_8uv.
\]

全P2はdelta_z、V_z及びPxyへquadratic(u,v)を加えるだけで式(1)を保つ。
Pxyはx ell(W)-y^2+y m(W)+q(W)であるが、以下の最弱証明にはその係数も不要。
normalizerで追加parameterを除いたchosen sliceではなく、complete K1と全九P2そのものを使う。

## 2. 一つの元projective Jacobi式

strict Jacobiを仮定せず、homogeneous representativeの元条件E wedge J=0を使う。
B内の全bracketが零なので、i,j in {z,u,v}に対し直接Leibnizから

\[
 J_{xij}=\delta_i A_j-\delta_j A_i,\qquad J_{zuv}=0.
\]

**Jxijは一般には零ではない。** {Bi,y}=-delta_iの寄与を必ず保持する。
元quadruple(x,z,u,v)のprojective Jacobiは、従って符号を除き

\[
 \det(E_B,A,\delta)=0,\qquad
 E_B=z\partial_z+u\partial_u+v\partial_v.\tag{2}
\]

これが全sourceの必要式である。
Euler補正を捨てたり、他のstrict Jacobiを余分に課したりしない。
Aの固有値又はnormal-form matrixのinvertibilityは不要。

## 3. actualDはx,yに依存しない

B内の全bracket零から、実Euler-Pfaffianは全体符号を除き

\[
 D_{\rm act}=
 \det\begin{pmatrix}
 z&u&v\\
 \{x,z\}&\{x,u\}&\{x,v\}\\
 \{y,z\}&\{y,u\}&\{y,v\}
 \end{pmatrix}
 =y\det(E_B,A,\delta)+\det(E_B,V,\delta).
\]

Pxyはこの式に現れない。式(2)によってy項が消え、残項はBのquintic。
xへの依存も初めから無い。
従って**全fullsolutionのactualDはC[z,u,v]だけに依存する**。
P2のparticular solutionは全九方向のまま含まれ、kernel torsorの取り違えはない。

irreducible quinticが二変数以下にしか依存することは複素数上不可能なので、指定部分ではessential数は正確に三。
既受理の `ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md`、magnetic-removal全source定理、
八curveの必要分類及び全ambient component/distinctness packetを適用する。
全pointは八23成分の閉和へ入り、G非優越も従う。
新しいsample curve、個別anchor又は高次obstructionを作る必要はない。

## 4. V1の不採用と完了範囲

旧未受理producer草案 `PLANE_M2_DOUBLE_TANGENT_CUBE_ORE_RECOGNITION.md`
（SHA 7ea7b8075ed523d77872a898abc069b1c3d151deb42fe0d6df3c2d9fcd6f9fbf）は、
Jxijを零とした段階が誤っていたため採用しない。旧bytesは履歴証拠として変更せず、
本V2を有効証明とする。結論は同じだが、正しい一Jacobi式によって直接示される。
18weightsによる旧論法は使わない。

本稿の読解が独立に受理されれば、有限jobで残ったh=L^3全源はE/Gに閉じる。
他三tangent模型の全像22以下、全空、whole moving像23以下と合わせ、
reduced tangent-conic m2 double全sourceのG非優越が揃う。
他三源の全E帰属、triangle/concurrent double全higher、nonreducedR、m3は別残余。
既code/report/job/20-file manifestは一切変更していない。
