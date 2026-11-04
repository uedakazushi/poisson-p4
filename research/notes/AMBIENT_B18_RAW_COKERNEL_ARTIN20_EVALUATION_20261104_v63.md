# B18 実 Artin20 source 評価と cokernel の閉ファイバー8

2026-11-04、v63。

実 conductor 積の閉ファイバー
\(A_{20}=R_0/(I_A I_B,\ell)\)、\(\ell=5(s_A-r_A)\) の正確な乗法表を使い、
両 normalizer の係数関数を有限代数上で評価した。実105単位方程式の全てを正確に満たす値、
実 physical 九座標の readback、両 \(\nu\)、十二組の \(a_A,a_B,b_A,b_B\) を保存した。

その値と実 raw26 商写像を合わせた結果、実比較 cokernel \(C\) について

\[
 \dim_{\mathbb Q} C/(\ell)=8
\]

を得た。312本の physical 作用列の rank は18である。
これは右完全な cokernel の特殊化であり、比較 kernel の特殊化ではない。
一般 DVR rank、kernel torsion、元の (D_0) の深さはこの結果から決めない。

## 1. 同じ actual source を有限代数で解く

入力の original source DAG は、CA pure48・CB pure48、実 physical 共通底九座標、
canonical cuts による両 actual normalizer である。
直前 preflight は、この両 pure48 が既受理 free4 normalizer packet のものと、
元 node ID の pruning map を通じて同一であることを確認している。

凍結済み joint105 の順序は CA57（pure48 と F9 の九式）、続いて CB48 である。
実 CA F9 の target は、正確に

\[
 (z_{A,0},\ldots,z_{A,6},r_A,r_A)\in A_{20}^{9}
\]

とする。raw physical 底を二次近似で代用しない。それら九座標も105未知数の一部として
有限代数の中で求め、元 physical rows
\(27,108,109,112,113,117,118,119,61\) の順に保存した。

未知数の中心は全て零。凍結原点 Jacobian の逆行列（866非零係数）を (J_0^{-1}\) とすると、

\[
 u\leftarrow u-J_0^{-1}(F(u)-\mathrm{target})
\]

により、残差の次数は少なくとも一つ上がる。\(\mathfrak m^5A_{20}=0\) なので、
零初期値から四回の補正を行い、その後の全105式が厳密に零になることを確認した。
各段階の105残差・補正・新未知数を保存した。これは新しい全九変数 formal jet solve ではなく、
確定した20次元有限代数における一回の単位方程式評価である。

必要な source 閉包18,584 nodes と六つの逆元を評価した。逆元はその中心が単位であることを
確認し、nilpotence による五項の幾何級数で求め、積が1であることも確認した。
元の定数 node は `['c', numerator, denominator]` として正確に評価している。
凍結 preflight の中心だけを読む箇所では分母を省略していたが、今回の正確な有理定数評価で
全式の中心と全六逆元の中心1を再確認した。source・表・順序の hash は証拠に固定している。

単位形式の一意性により、この105解は actual 正則 normalizer 関数の (A_{20}\) への制限である。
古い非単位 (r/g\) 商を新たに評価したり、その未知微分を零と置いたりしていない。

## 2. 実 raw26 への作用

別の凍結 packet `AMBIENT_B18_RAW_COKER26_CLOSED_ACTIONS_20261104_v63` は

\[
 C_{\mathrm{raw}}/(\ell)=
 (A_{20}/I_A)\tau_A\oplus(A_{20}/I_B)\tau_B
 \oplus A_{20}\tau_A\tau_B
\]

の実26基底を与える。最初の両三次元 sector の基底は共に CA の \((1,r_A,r_A^2)\) である。
残りは凍結20基底である。二つの20→3商写像と共通代表持上げをそのまま利用した。

十二個の actual physical differences は

\[
 F_j=\delta_j+b_{A,j}\tau_A-b_{B,j}\tau_B,\qquad
 \delta_j=a_{A,j}-a_{B,j}
\]

であり、source の元 node ID と全 (A_{20}\) 値を保存した。
各作用は商写像・代表持上げを伴う block 表示

\[
 \begin{pmatrix}
 \delta_j&0&-\nu_Bb_{B,j}\\
 0&\delta_j&\nu_Ab_{A,j}\\
 -b_{B,j}&b_{A,j}&\delta_j
 \end{pmatrix}
\]

から構成した。
両 actual conductor 八生成元の商への消滅、全64 actual source 積の (A_{20}\) における消滅、
\(I_A b_B=I_B b_A=0\) を確認した。従って代表持上げの選択によらず作用が定まる。
十二行列は八つの座標作用と可換であり、相互にも可換であることを全係数で確認した。

## 3. 実 cokernel の右完全な特殊化

既受理の actual tensor/cokernel 表示は

\[
 C=C_{\rm raw}/\sum_{j=1}^{12}F_j C_{\rm raw}
\]

である。従って右完全性だけから、閉ファイバーも同じ十二作用の像の和で割った商となる。
全312列について有理数上の pivot 証拠を保存し、その rank が18、商が8次元であることを確認した。
実 raw26 基底の位置 \(0,1,2,3,4,5,6,14\) が商基底となる。
八座標の誘導作用行列と26→8商写像も保存した。

この計算は、特殊化行列の kernel が比較 kernel の特殊化であると仮定していない。
有限 DVR module としての一般 rank・torsion、必要な完全比較 complex は引き続き別の入力である。

## 4. 証拠と実行記録

主要証拠は同名 report directory の次の二ファイルである。

- `actual-source105-Artin20-values.json`：全 source 値、四補正、全105式、六逆元、physical readback。
- `actual-physical12-actions-and-C-special.json`：十二作用、全 intrinsic checks、312列と pivot、8次元商作用。

source 評価は0.54秒。最初の実行では action の計算後、相対 CLI path を絶対 path の
`relative_to` へ渡した metadata 保存で停止した。その source・失敗 prefix を保存し、
path の正規化だけを修正した。保存済み source 値から action 保存だけを0.037秒で完了し、
105 solve は再実行していない。最大 RSS は約24 MB、出力は全体で十分32 MB以内である。

初回用 source は `worker-at-first-launch.py`、失敗 metadata は
`first-launch-metadata-error.json` に保存した。通常の再現入口は

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_raw_cokernel_Artin20_evaluation_20261104_v63.py --raw26-evidence reports/AMBIENT_B18_RAW_COKER26_CLOSED_ACTIONS_20261104_v63/actual-raw26-actions-and-quotients.json
```

であり、修正後は一回で両段階を保存する。
