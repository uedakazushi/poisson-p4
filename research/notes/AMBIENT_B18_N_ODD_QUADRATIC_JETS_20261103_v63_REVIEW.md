# generic B18のactual N odd二次項：独立検証

`AMBIENT_B18_N_ODD_QUADRATIC_JETS_20261103_v63` の有限source計算と、actual odd idealのminimal生成数9以上という帰結を受理する。全N像の正確なidealへのNakayamaによる応用は別のproof companionとして区別する。

## 全sourceとgauge

raw CA式を直接canonical化し、必要な座標交換を行い、固定旗N19基底を別実装で再構成した。全16個のparabolic生成子が全19基底を保つこと、中心でのparabolic orbitの10 minor、残る9個のfiber方向の補空間を照合した。全9個のflag方向の旗微分のrank9と、fiber19に旗を加えたaugmented sourceのrank28も確認した。したがって、rank2点だけに使える五つの旗変数で全sourceを代用していない。

実際の二変数planeは

\[
P(W,V)=P_0+W B_2+V(5B_1/2+11B_8/2+B_{16})
\]

である。このplane上の全150×28 source Jacobianを、19固定fiber列と9物理GL Lie列から独立に作り直した。全entryが保存証明書と一致する。parabolic10とfiber補空間9のétale gaugeは全sourceに関する既受理のcritical idealを適用するためにも確認した。

## actual baseを保つ微分の二次まで

列22以外をB、列22をbとし、選択した27 actual target行をrとする。\(B_r=T_0+T_1\) に対し、別の実装では直接

\[
R_0=T_0^{-1},\quad R_1=-R_0T_1R_0,\quad
R_2=R_0T_1R_0T_1R_0
\]

を構成した。\((R_0+R_1+R_2)(T_0+T_1)=1\) が次数2まで成立することを全係数で確認し、\(-R b\) の次数0、1、2を抽出した。producerの再帰式とは別にこのinverse展開から計算した三つのvectorは、保存された全27係数と一致する。

その後、全150 target微分の次数0、1が零、次数2が保存vectorに一致することを確認した。選択27行は二次まで零であり、残るbranch row27も二次まで零である。

## Schur二次項をactual odd二次項として扱える理由

共通actual底は \(R=\mathbf C[[Y_{28}]]\)、全normalizerは \(B=R[t]/(t^2-\nu)\) である。上のvector field Dは27 base座標を保ち、kernel flag coordinateの成分を1に正規化しているので、\(Dt(0)=1\) である。

branch row27の微分は、\(\nu=t^2\) と \(\partial\nu/\partial Y_{27}\ne0\) から、単位倍の \(2tDt\) である。その次数1、2がともに零なので、plane上でtは次数3以上になる。各target出力 \(F_i=a_i(R)+b_i(R)t\) にDを作用させると、二次以下ではeven部分やtを伴う補正が消え、\(DF_i=Dt\,b_i\) となる。

planeのRへの一次微分はrank2で、全actual odd一次行列Aを掛けると零である。従って各b_iのplaneへの一次制限も零であり、\(Dt(0)=1\) だけでSchur二次項はactual b_iの二次制限と一致する。古いcore pと今回のflag pを同じ座標と仮定する必要はない。同じrank2 R-algebraの二つの基底1,pはR上のaffine linear変換で関係し、tracefree tは中心で値1の単位倍で変わる。full R firstjetが零のodd生成元の二次classは、この変更でも不変である。

## 七つの一次classと二つの追加class

既に独立検証したactual Aの全150×28 entryを固定し、今回そのrank7と基底行

\[
18,19,28,29,60,61,67
\]

を再計算した。全行の一次結合と \(d\nu=-db_{67}/8\) も照合した。

全b_iからこれら七行の定数一次結合を引くと、**full R上の一次項が零**の元を得る。独立に再構成した二次制限は

\[
g_{46}=-4WV,\qquad g_{72}=-6V^2
\]

であり、係数minorは24である。全150 residual二次式も一致した。任意の \(\mathfrak m_R I\) の二次項は、一次idealのkernelに含まれるこのplaneへ制限すれば零になる。従って二つのclassは \(I/\mathfrak m_R I\) で独立で、七つの一次classに加えて \(\mu(I)\ge9\) を与える。

選んだ \(-b_{67}/8\) のplane二次項は \(-9V^2/4\) である。ただしzero-firstjet生成元を加えれば変わるため、これをν自体のcanonicalな二次項とは呼ばない。

この結果はpaired source arcの不存在から生成数を推測するものではなく、actual targetのodd係数そのものの計算である。従って二つの八生成ansatzはどちらも排除するが、任意の高接触次数のpaired arcの不存在を証明したものではない。元のfull9旗のpaired三次計算は、この結論の依存関係に含めない。

検証codeは `code/ambient_B18_N_odd2_review_20261103_v63.py`、receiptは `reports/AMBIENT_B18_N_ODD_QUADRATIC_JETS_20261103_v63_REVIEW.json` である。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/ambient_B18_N_odd2_review_20261103_v63.py \
  --replay work/B18-N-odd2-producer-replay \
  --output work/B18-N-odd2-independent-review.json
```
