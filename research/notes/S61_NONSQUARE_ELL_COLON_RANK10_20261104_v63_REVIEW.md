# ell-colon / cubic11 と保存点 rank10 の独立監査

2026-11-04。対象は凍結済み `S61_NONSQUARE_ELL_COLON_CUBIC11_20261104_v63`
および `S61_NONSQUARE_ELL_CUBIC_RANK10_20261104_v63`。親の frozen source を変更せず、
保存済み有理点 `(a,e,d)=(1,2,3)` のみを使う。独立実行の数値・hash・全係数の一致は
本 note と同名の REVIEW report/certificate に記録する。新点、素数、GB、関数体 solver、
四次式全展開は使わない。

## 1. degree3 colon と 16/64/48

FULLU の base ring を R、七 light 変数の多項式環を S とし、
J0=(q0,q1,q3,q5)、I=J0+(ell Z,ell W) と書く。
core の既知因子 J と ideal J0 を区別する。ell の x1 係数は既知 unit なので、
ell は S の多項式座標にできる。

U78 は mod ell の四 quad の degree2-cofactor map 84→126 が rank78 の開である。
既存の四 quad coefficient dual により六 Koszul directions は split subbundle である。
rank78 の開では kernel も rank6 bundle であり、各 fiber 上で一致するから、この六方向が
**whole kernel module** になる。単なる generic kernel equality ではない。

ell f3=Σ qi Ai2 とする。mod ell の係数四組は上の六 Koszul の組合せなので、その組合せを
S 上の六 Koszul で持ち上げて引けば、各 Ai が ell で割れる。ell を cancel するのは
多項式環 S 上であり、S/J0 上での非零因子性を仮定していない。従って

    (J0:ell)3=(J0)3,
    (I:ell)3=(J0)3+(Z,W)S2.

標準的な homogeneous colon exact sequence の degree4 部分は

    0 → V16=(S/(J0,Z,W))3 --ell→ B64=(S/I)4
                                      → V48=(S/(I,ell))4 → 0.

左の16は、受理済みの実際の五変数 product map 20→35 の rank19 による **35−19**。
四 quad CI からの15ではない。右は126−78=48。左右が constant-rank quotient bundle
なので中央も rank64 bundle になる。全次数の ell regularity、CI、unmixedness、CM は使わない。
この議論は U78 上であり、その外の有限点をこの開へ勝手に加えない。

## 2. actual C122 内の canonical11 と正確な lift

heavy7 chart で実 genuine16 cubic を、七式 fi=H xi−pi と九 pure 式 rj に可逆に変換する。
Z=3v0−v1+v6、W=v3 とし、fZ=3fv0−fv1+fv6、fW=fv3 とする。
actual quartic source 上の十一の lift は

    ell*rj (j=1,…,9),
    KZ=H*(ell Z)−ell*fZ=ell*pZ,
    KW=H*(ell W)−ell*fW=ell*pW.

後二つは Hq source を重い正規形へ送る際の固定 commutator 補正を含む。
単独の NF_H(H ell Z)=ell pZ を literal identity として使わない。
受理済み 122 recipe と aligned source row/frame の可逆輸送により同じ式を C122 へ送る。
十五の full-kernel directions や三十二 genuine quartic residual は保持する。

九 ell*r source は九つの別個の pure-cubic/light block に属し、ell_x1 unit で split する。
残り二つは Hq block の ell Z,ell W にあり、六 quad の prime-basis row change det−80 により
split 独立である。pure-quad Schur と35行削除はこの二つの Hq component を失わせない。
従って十一 **source** は split 独立である。その十一の **target image** の独立性とは別である。

Genuine16 を Z=W=0 に制限すると H-light target は五次元。
元の七 heavy directions のうち Z,W の二つが消え、九 pure directions とともに kernel11 を与える。
quad products の十九方向を35 pure cubicsから割ると、上の canonical11 の V16 への像と
同じ intrinsic map 11→16 になる。可逆 source/target chart change の下の一致であり、
保存点で一致したことだけからパラメータ全体の同一性を推論していない。

## 3. 独立係数・Schur・minor 照合

保存された57×91 actual cubic matrix と元 quadratic sparse coefficientsから再構成する。
author の binomial restriction helper を使わず、変数像の反復 sparse convolution で
v3=0、v6=x1−3x0 を代入し、ell 自体は自由変数として残す。
四 quad×五変数の20×35、genuine16 の16×40、元 heavy16×7 の **全1452有理係数** を
native input と照合する。source monomial order、元 row index、既存 polynomial source hash
および parametric circuit の保存 pivot index を検証する。

author の二段階 Schur を再利用せず、quad19 と heavy5 の24×24 block を一括逆行列にして、
残る11×16の176係数をすべて独立に計算する。二十 quad 行の Schur residual は全て零。
既に選ばれた19、5、7、10、34次 minor の値は、各行の分母を払った整数 Bareiss 法で照合する。
別 minor の探索はしない。11行が保存済み十行から正確に再構成されることも確認する。

結果はこの保存点で rank **10**。parametric 34次 minor は同じ実 polynomial source の minor で、
その点での値は det19*det5*det10≠0 である。従って generic rank≥10 は証明される。
rank≤10 は **この点だけ**。generic 上界は source rank の11のままであり、[10,11]を保つ。
Delta10、quad19、heavy5 の minor は追加 chart factor で、FULLU の既知 unit へ加えない。

## 4. 条件付き E49→V6 の同じ cokernel

U78・heavy7・Delta10 chart 上で十一のうち十の image directions を split する。
十五の既知 kernel directions とは target image が零/独立なので交わり零。
source122 を十五 kernel と十 image source で割ると rank97、target64 を十 image で割ると rank54。
target は exact extension 0→V6→B54→V48→0 を持つ。
さらに **actual source97→V48 が全射** の開では、その kernel E49 と V6 の map が
元の full quartic cokernel を保つ。これは diagram/snake lemma による bundle criterion であり、
固定49×6係数行列を既に展開したとの主張ではない。全射条件は、この REVIEW で証明しない。
十方向消去前の同じ議論は E59→V16、十一 image が独立なら E48→V5 を与えるが、
generic rank11 を今回の保存点から主張しない。

したがって受理範囲は degree3 colon、正確な canonical11 source recipe、保存点 exact rank10、
generic rank≥10 と上の条件付き同 cokernel criterion。全点57、全四次239、全点 unit、
generic upper10、Delta10 上だけでの projected V48 全射は未証明のまま保持する。

再実行（新しい空 output directory を指定）:

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_ell_colon_rank10_review_20261104_v63.py --output work/marked-square-v63/ell-colon-rank10-review-REPLAY
```
