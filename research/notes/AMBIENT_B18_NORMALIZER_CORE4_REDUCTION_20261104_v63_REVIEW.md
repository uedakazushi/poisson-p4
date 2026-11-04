# Actual E4 core4 typed reductionの独立監査

2026-11-04。別の独立N-pairレビューを受理gateとして使う。
author runnerを実行せず、元pure48 source nodeから全98式・十二physical差を独立に再構築した。
元source DAG41315ノードのprefix一致と全suffix演算を確認し、CA/CB各48 unit未知数と
rho27,rho61のeffective Schur minor -20/9を照合した。

N_BETAは、既受理のnormal source上の一意なregular trace rootsとrank8 inverseから定義される。
中心でconductorをunitとして反転していない。
L=6*(tauB-N_BETA)の一次項を独立N-pair normalと照合し、全十二FのL-regular除法を受理した。

cached二次式をこのL一次項で**独自に係数除法**し、全45 monomialへ掛け戻して全十二linear quotientを検証した。
五行18,19,60,28,29と五座標108,109,113,118,119のmatrixは

```
diag(1/12,1/12,1/24,1/12,1/12), det=1/497664.
```

core4=(x112,x117,tauA,tauB)。七retained rowsは67,72,63,64,68,73,74である。
SCHUR98、N_BETA、L、G12、CORE5、K7の全六typed operator interfaceを元sourceに結んだ。
これは高次係数展開や、全operatorを多項式unit系だけに置き換えた証明ではない。

S9でJ=(G12)と書けば、formal座標を五つのselected quotientとcore4に変えることで
J=(U5,K7(core4))となる。元E4 idealはL'*Jのままであり、五つのL'*U_i関係を残す。
変換後のdL'/d tauBは原点で6、従ってunitである。
その全関数が常に定数6だとは解釈しない。このclarificationはauthorにも確認済み。

S9はdomain、Lは非零だから

```
0 -> S9/J --L--> S9/(LJ) -> S9/(L) -> 0
```

は完全。S9/(L)はregular8なので、H0(E4)=H0(S9/J)=H0(core4/K7)を得る。
E4そのものをcore4 quotientと同一視しない。
またD0→E4の比較kernelは未処理なので、D0のsocleやdepthへこの等号を移さない。

独立code `code/ambient_B18_normalizer_core4_review_20261104_v63.py`。
Report/folder `AMBIENT_B18_NORMALIZER_CORE4_REDUCTION_20261104_v63_REVIEW`。
新source solve、高次jet、GB、new point、消える分母のunit inversionは無い。
