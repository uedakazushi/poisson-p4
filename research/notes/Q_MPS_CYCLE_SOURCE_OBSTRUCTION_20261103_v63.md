# T39/T40に対する実extended toric arc源の厳密な制約

2026-11-03。保存済みextended arcが、相異なる固有値とrho(rho²−1)≠0の開で
T39またはT40の純toric periodへ入るなら、どの120標識置換でもそれぞれ
**a=2/9、a=−8/9**に限られる。
これは特定の実arc源への必要条件であり、対象成分全体への非所属を示さない。
特にa=0についての非所属は一切結論しない。

## 1. Vandermondeの補小行列による低次数化

d=(0,1,z,x,y)、Vの(i,k)成分をd_i^k (i,k=1,…,4)とする。
Delta=det(V)=zxy(z−1)(x−1)(y−1)(x−z)(y−z)(y−x)。
標識置換した対象periodを線形Pfaffian双対でbivector行列Aへ戻す。
H=Delta*V^(-1)*A*V^(-T)の六成分はVの2×2補小行列の線形和である。
従って巨大な有理逆行列を作らず、三変数の小さい多項式で計算できる。
producerは全六成分についてV H V^T=Delta AをQQ上で確認する。

extended arcのpower basis係数は(alpha,beta,rho,1,0,0)だから、
対象の非零period lineへ入るにはH24=H34=0、H23≠0が必要であり、
rho=H14/H23となる。

## 2. 全標識に対する理想所属恒等式

各対象の120標識は12個のprojective period formに分かれる。
その全24場合について、実際のQQ三変数多項式U,Vを保存し、

    9*a0*(H14²−H23²)−8*H14*H23=U*H24+V*H34

を厳密に再乗算して検査した。
a0はT39で2/9、T40で−8/9である。
従って当該source上では9a0(rho²−1)=8rhoとなり、
a=8rho/[9(rho²−1)]=a0を得る。
単なるGroebner normal form零や有限体標本ではなく、元の二方程式に対する明示cofactor証明である。

## 3. 結論の範囲と保存

T39/T40への既知の実所属a=2/9,−8/9は、歴史的なS61 quadratic graphの別の証明による。
このpacketはその所属を新たに導くものではなく、今回のextended arc sourceが
一般rhoを動かせない理由を正確に示す。
他の退化族、smoothing方向、a=0のrank2境界、対象成分全体のQa逆像はこの必要条件の射程外である。

producer `code/q_mps_cycle_source_obstruction_20261103_v63.py`、
cofactor `code/data/q_mps_cycle_source_cofactors_20261103_v63.json`、
receipt `reports/Q_MPS_CYCLE_SOURCE_OBSTRUCTION_20261103_v63/`。
source制約の保存後は、一般S61の未知支持の研究へ戻る。
