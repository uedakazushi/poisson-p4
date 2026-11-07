# tangent normal-double・h=L³の全源：一projective Jacobi式による別作者認識

2026-11-07、`/root/quartic_review`。**全十八parameter sourceの独立恒等式確認と全支持の数学はPASS。**
actualD被約既約quinticの指定rank四部分は、全点が既essential-three quintic八成分へ入る。
既四源reviewと合わせ、reduced tangent-conic・m二double-leadingの全source像は二十三以下となる。
他三源の全E帰属、triangle/concurrent double全higher又は全M2分類を主張しない。

## 1. 採用する証明と不採用の旧草案

有効producer証明は `PLANE_M2_DOUBLE_TANGENT_CUBE_ORE_RECOGNITION_V2.md`
（fd6e7bf5…）である。
旧7ea7b807…草案はstrict JxBiBjを零とした段階が誤りだったため**採用しない**。
実際JxBiBjには{Bi,y}の寄与が残る。
旧bytesは変更せず、18個のad-weightや誤ったradial inferenceを依存に含めない。

全sourceの入口は不変の `CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_20261107.json`
（ebcafe36…）constant_hcubeと、既独立四源review（966a851d…）。
全350×30の完全K1九と全P2九を使う。
**P2をkerJ2三に狭めない。** 非同次weight二式のparticular solutionも落とさず、全九方向を許す。
以前の粗27というreport及び四源partial reviewはそのscopeのまま不変に保つ。

## 2. 全十八係数のsource形

x,y,z,u,vの順に座標を置き、B=C[z,u,v]とする。
全K1をn0,...,n8、全P2をr0,...,r8で戻すと、

`{z,u}={z,v}={u,v}=0`、
`{y,B}=delta(B)`、`{x,B}=y A(B)+V(B)`

である。Aはlinear、delta,Vはquadraticで、全てx,y非依存。
complete kernelと全九P2から直接組み立て、各monomialを確認した。

具体的に

`A=(2z+(9n6+3n8)u+(2n2+4n4)v, −16u, −10v)`、

`delta_z=z²+n2*z*v+(9/2)(n6+n8)z*u+r6*u²+r7*u*v+r8*v²`、
`delta_u=−8zu+n4uv`、`delta_v=−5zv+n8uv`。

P2のxz成分はV_zへ、yz成分はdelta_zへ、xy成分はPxyへ入るだけなので、
全P2九を戻しても上のsource形を保つ。
Pxyの係数又は−y²項に余分な制約を課す必要もない。
additional normalizerで選んだsliceではなく、全必要source十八係数を同時に扱っている。

## 3. strict式を勝手に零にしない直接証明

W=(z,u,v)、E_B=(z,u,v)と書く。
B内のbracket零と{y,B}=delta、{x,B}=yA+Vから、i,j in Bには

`J_xij=delta_i A_j−delta_j A_i`、`J_zuv=0`

がLeibnizで従う。
ここでJはsingle cyclic strict Jacobiであり、J_xijは一般には零でない。
元のprojective JacobiはE wedge J零だけを要求する。

quadruple **(x,z,u,v)** の一式を展開すると、正確に

`Jac_projective_(x,z,u,v)=det(E_B,A,delta)`。

一方Euler-signed Pfaffianを全五項から計算すると、B内のbracket零から

`D_actual=−det(E_B,{x,B},delta)`
`        =−y det(E_B,A,delta)−det(E_B,V,delta)`。

PxyはこのactualD式へ寄与しない。
従って全十八parameter上で、条件を課す前から

**`partial_x D_actual=0`**、
**`partial_y D_actual=−Jac_projective_(x,z,u,v)`**

がpolynomial identityとして成立する。
この正しいprojective一式だけで、全Poisson solutionのD_actualはC[z,u,v]に入る。
他のstrict Jacobiを零と置いたり、Aの固有値・invertibility・generic diagonalizabilityを仮定したりしない。
source上の全period、zero coefficient、全higher compatibilityを残す。

## 4. 独立有限確認と最弱の全E/G含意

新checker `classification_plane_m2_double_cube_review.py` はproducer/CASをimportせず、
既独立のFraction sparse演算で上の全K1九＋P2九を再構成した。
有限job `jobs/plane_m2_double_tangent_cube_separate_review.json` gate PASS後、
**0.049秒、25,036KiB**で全source構造、全A/delta係数、元Jacobi一式、actualPf全五項と両微分恒等式を確認した。
新report `CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_CUBE_SEPARATE_REVIEW_20261107.json`
に実行reviewer全文と全input hashを保持する。
以前の四源／全七壁計算は再実行しない。

複素数上の既約quinticがessential二以下になることはbinary factorizationに反するので、
指定rank四・既約部分のessential数は正確に三。
既独立受理済みのessential-three Saito全source・全磁気除去・八curve完全分類・
全八ambient componenthood/distinctness定理を適用出来る。
従ってこのh=L³全源の指定部分は八既23成分の閉和に入り、G像も二十三以下。
可約又は非被約actualDの全E帰属をここから追加しない。

前の四源partial reviewのhf/Z零像≤22、hf/ZN全空、moving全t非零像≤23と合わせると、
**tangent-conic reduced R・m二normal-double全枝のG像は二十三以下**。
先の粗27は不正な数値だったのではなく、今回の構造恒等式より弱い上界だった。
全normal-double又は全G/Eを完了したとは呼ばない。
