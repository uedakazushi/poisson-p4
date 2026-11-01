# 自由 M0 旗端点：実際の R-adic 有限商と加群準備の WIP

2026年11月2日。追加の高次 native 計算を止め、既存の成果と失敗条件を保存する。
本ノートは一つの有限体点での準備計算であり、自由 M0 層の飽和除去を主張しない。

## 1. 点、filtration、得られた有限商

対象は既存 FLAG_ENDPOINT_TABLE の free_M0 点である。
通常の16個の一次座標では

    (0,0,0,0,−1/3,6,4,7,8,5,9,10,11,1/3,−23,17)

で、旗 Q=−x²∂y∧∂w の chart u38=−1 を保つ。
INITIAL_PROBE の11接座標 v0,…,v10 を用い、
normal 変数を v0,…,v7、base を u=(v8,v9,v10) とする。
これは universal normal-form の a=(1,2,3,4,5) 点とは別である。
係数体は F_32003。S は元135式の指定点での完備局所環、a=(u)S と書く。
ここで扱う filtration は最大 ideal filtration ではなく a-adic filtration である。

元135式の陰関数展開に base の冪を加えた三つの有限商で、保存された計算は

| 商 | 長さ | 最大 ideal の零冪 |
| --- | ---: | ---: |
| S/a | 141 | m^8=0 |
| S/a² | 542 | m^9=0 |
| S/a³ | 1303 | m^10=0 |

を与えた。最初の商の Hilbert 関数は (1,8,24,38,38,23,8,1)。
それぞれ次数10、8、9までの元式から計算している。
保存 initial leading ideal が示した零冪は切断境界だけからの推測ではない。
その次数の全単項式が、境界より低次の実在初期式の leading monomial で割り切れる。
従って実際の形式環でも m^e⊂m^(e+1) となり、Nakayama により m^e=0。
このため、当該有限体の Artin 商については有限切断で完結している。
ただし raw135 展開・native standard basis 自体は本 snapshot validator では再実行しない。

## 2. 141生成への圧縮

横断商の141標準単項式は complete Nakayama により R=F_32003[[u0,u1,u2]] 上の生成候補となる。
保存された正常形から、base 一次423列の rank は401、kernel は22次元。
base 二次846列の rank は761、kernel は85次元である。
一次関係の変数倍は60次元なので、新しい二次関係は25本。
よって実際の R-adic 最初の三層の寸法は

    141,401,761.

22本の一次関係だけの P=F_32003[u0,u1,u2] 加群の rank は125、
85本の二次関係を加えたものは108だった。後者の初期 Hilbert 値は
(141,401,761,1229,1805,2489)。これらの関数体 rank は保存 M2 実行記録であり、
snapshot validator は M2 を再実行せず、数値と入力を保存する。
validator が新たに点検するのは、全正常形行列の rank401/761、
22/85個の kernel ベクトル、独立な一次倍60個である。

同じ実際の S/a³ 正常形で t は非零（保存正常形1053項）だが t² は零だった。
従って得た情報はこの有限体点での t²∈a³S だけである。
全 S、標数零、relative β のいずれについても t²=0 は導いていない。
この t² の native 還元は計算記録として保存し、今回の小さい validator で別実装したとはしない。

## 3. 全局所環へ進むために必要な有限橋

archive 中の R_ADIC_TORSION_FREE_RECOGNITION 草稿は条件付き補題である。
有限 R 代数 S と実在閉商 S→D があり、P 上の有限表示 M が

    M → gr_a S → gr_a D

へ全射するとする。rank_P M=rank_R S、
H=ker(M→gr_a D) が P-torsion-free、t^h=0 が D と S⊗Frac(R) の双方で成り立つなら、
t^h=0 が S で従う。
実際 gr_a S の P-rank は R 上の重複度、従って rank_R S に等しく、
M→gr_a S の kernel は rank 零で H の部分加群なので零。
S→D は a-adic filtration に strict である。
N=ker(S→D) に誘導した filtration の gr が H となるため、separatedness から N は R-torsion-free。
t^h は N の torsion 元なので消える。

現状の候補 M は rank108であり、期待される既知枝の総 rank70に達していない。
H の torsion-freeness も未証明である。
既知枝の下界をこの有限体に移す仮定、QQ への移行、relative な飽和支持の議論も別途必要。
従ってこの補題の適用を完了したという主張はしない。
旧 m-adic initial 商の dim2/degree31 の差を、実際の R-torsion と同定もしない。

次の改善候補として、既知 homogeneous lower initial ideal J を P 上の有限加群にし、
G=F[v]/J の G/(u) が最高次数 e−1 を持つことを確認する方法がある。
その場合 graded Nakayama から全 r に対し m^(r+e−1)⊂a^r が従う。
こちらでは J10/(u) の有限性と最高次数をまだ確認していない。
actual 横断141だけでこの前提を代用しない。

## 4. 再現資材と停止した計算

小型 archive は

    research/wip/S61_MINIMAL_FREE_FLAG_R_ADIC_PREPARATION_20261102.tar.gz

で、元の探索ソース、ログ、疎 initial 基底、141/564/1410列の正常形情報、
加群行列、再構成順序の README と SHA256 manifest を含む。
320MB の raw 入力、74MB の完全 filtered 基底、native バイナリ、ライブラリは含めない。
入力は追跡済みの135式と endpoint report から再生成できる。
探索ソースは元の /tmp パスと M2 コマンドを保持しており、統一された production runner ではない。

再利用できる二つの native helper も保存する。
`truncated_local_module_basis_20261102.cpp` は既存の local standard basis に query 還元と完全基底出力を加え、
`truncated_local_module_normal_form_20261102.cpp` は保存基底による query 正常形だけを計算する。
後者の1410列は元の出力と byte 単位で一致した。
今回それらを改めて compile/run することはしない。

    .venv/bin/python code/s61_minimal_free_flag_r_adic_preparation_20261102.py --output reports/S61_MINIMAL_FREE_FLAG_R_ADIC_PREPARATION_20261102.json

このコマンドは archive の integrity、leading staircases と小行列の整合を検証する。
原135式からの全再計算と混同しない。

S/a⁴ を狙った raw10 展開は138.45秒で作れたが、standard basis は元順序240秒で未完。
base の15個の四次単項式を先に置いた異なる順序も180秒で未完だった。
base² を明示的に加えた local Mora の別試行は90秒で未完。
これらを同条件で延長・反復しない。
既知の有限商での bitset、flat hash、parameter quotient mask の最適化は大幅な改善を生まず、
次手は141生成の加群を直接扱う方法である。

全 initial ideal 等号、actual torsion の認識、ambient t 冪零、自由 M0 層の除去、
一般 affine fiber の有限上界、S61 の成分尽くしは、いずれも本保存では未解決である。
