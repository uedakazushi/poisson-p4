# actual free281 有限 circuit の独立監査

2026-11-04。対象は凍結済み `S61_NONSQUARE_FREE281_CIRCUIT_20261104_v63`。
base は FULLU=D(adFBCJ)。直前に受理した compact kernel38 と既存 P22 literal identity
を再利用するが、author script は再実行しない。source319/304/282 の全 index、
residual4とP22の composition、最後の既知 unit 消去を独立な expression proof checker で照合する。

## 1. source304 への正確な輸送

既存 atlas の retained304 と deleted15 が319 source indexをちょうど分割し、
先頭174が元の六 quad×29 monomial であることを確認する。
compact vector N は quartic33 と residual4の先の三 pivot座標176,177,237で零。
従って source319→304 は

    N' = N − N[181] Rlast

だけでよい。Rlast は、受理済み D15 の actual319 roots に residual-transform U と
正規化 pivot-transform T の最終行を合成したもの。
quadratic174 componentについてもこの合成を再構成する。
cubic componentは受理済み post-pivot matrix、quartic componentは受理済み exact zero を使う。
全319 Rlast roots、全319 N' rootsが保存 circuit と一致し、deleted15座標がすべて零になることを確認する。
N[181] は単に掛ける coefficient であり、これを反転しない。

## 2. Bnew/Cnew と、任意の B22 に対する P22 complement

旧 arbitrary B22 は separatorの元 source319座標182,186を両方含む。
従って旧282座標でこれらがそのまま残るとはいえない。
新しい B22 は retained304 surplus130から182,186を除いた最初の22本。
残り C108には182,186が残り、新282の座標174,175になる。
全 partition、順序、source319の定義との一致を検証する。

source304を (x174,u22,h108) とし、元 quadratic sourceの22 kernel columnsを K、
受理済み left inverseを P と書く。AK=0、PK=I22 は元 source 上の literal identity。
Bnewの選択に依存する仮定はない。

    k=Px−u,  z=x−Kk,  h=h

の逆変換は

    x=Kk+z,  u=Pz,  h=h.

実際、Pz=Px−PK(Px−u)=uである。よって304 sourceの自由な座標変換であり、
source map M304=[A,Bnew,Cnew] は

    M304(x,u,h)=AKk+(A+Bnew P)z+Cnew h
              =(A+Bnew P)z+Cnew h

となる。このため任意の Bnewで自由282表示 [A+Bnew P,Cnew] が得られる。
Kの22方向を source surplus の−Iと取り違えず、この全 source automorphism を使う。

レビューでは pinned P22 circuit の174 input variablesと線形 node topologyを検証し、
N' の x174をその全 output22へ代入する。残りの k22 と Kk174も元 sparseKから構成する。
P22 composition、k22、z174+h108の全282 roots が保存 circuit と構造的に一致する。
これは全 root を巨大な core 多項式へ展開せずに済む exact composition proofであり、点での比較ではない。

## 3. 唯一の scalar equality と unit coordinate 消去

新282 kernel の p=174,q=175は literal Cnew coordinates。
compact certificate の N182+N186=theta と Rlast182+Rlast186=0から

    v_p+v_q=theta=−1024000*a*C³*J/243

が従う。レビュー native は保存 scalar DAGのこの identityを sparse QQ numerator/denominator arithmeticで
独立に直接計算し、零 numeratorを確認する。他の大きな root は展開しない。
全898 atomの denominatorについても、非零有理数と既知 a,d,F,B,C,Jだけを順に厳密に取り除き、
未知 factor が残らないことを確認する。全10236 nodeの topological referenceも確認する。
theta自体が上の既知 unit polynomialであることを別に照合する。

係数座標の定数変換 w_p=v_p+v_q,w_q=v_q はdet=1。
generatorは M'_p=M_p、M'_q=M_q−M_pとなる。
変換後 kernel の p成分はthetaなので、thetaで正規化した kernelと e_i(i≠p) は
自由282 sourceの基底をなす。p方向は source map で零になり、残り281 generatorの像は
元の像と等しい。従って **全 FULLU 上で自由281→自由239の有限 recipeが同じ actual cokernelを与える**。

281 sourceの元319への liftも正しい。まず w_p=0を挿入し、v_p=−w_q,v_q=w_qへ戻し、
v=(z174,h108) として304 sourceへ (z,Pz,h) を入れ、deleted319 coordinatesを零で埋める。
これは上の generator recipe と同じ線形 map である。

## 4. scope と実行

新しい kernel方向の数は38のままで、目標80まで残り42。
free281は固定した元281行の部分行列ではなく、保存された source linear combinations による表示。
281×239の全係数を展開したとの主張ではない。全四次239、全点 cubic57、全 kernel exhaustion、
target cokernel零は未証明。本監査で新点・素数・rank・minor・GBは使わない。

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_free281_circuit_review_20261104_v63.py --output work/marked-square-v63/free281-review-REPLAY
```
