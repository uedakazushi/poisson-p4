# 全評価fiberとS61によるT39/T40所属判定の独立確認

2026-11-07、PRIMARY E。rootはproducerの全幾何と原S61/Q51/Q52のscopeを読解し受理した。
同時に、既E6 smooth-flexの単一faithful chartを異なるFraction実装で全係数確認した。
旧診断のbytesは不変であり、その否定の射程を新しい幾何定理で強める。

一般rank四の任意の全projective既約成分Cを取り、そのaffineconeをXとする。
connected PGL5は有限個の成分それぞれを保つ。pでの非退化値openの逆像X_Uは非空既約。
GL4は六次元の非退化skew値Uへ推移的、Jのstabilizer Sp4はconnectedなので、
X_U=GL4×^Sp4Y、Y=X∩{ev_p=J}の全fiberは既約である。
実際、Yの各既約成分をconnected Sp4が保ち、二つあればassociatedbundleにも二つの閉成分が降りる。
これはgeneric smoothfiberだけの主張ではなく全fiberの被約台の議論である。
dimC=dならdimY=d−5となる。

pを固定し接写像恒等の実perspective Ga4、z↦z/(1+az)は全Jfiberを保つ。
一定Jのpushforwardの一次項は
T_a=−2(az)J+z(Ja)^t−(Ja)z^tである。
保存Euler gaugeのP0iは定数/一次係数が零なので、四raw pivots61..64はintrinsic L12の四係数に等しい。
T12=−3a1z1−3a2z2−2a3z3−2a4z4、response det36は定数単元。
従ってa=−diag(−3,−3,−2,−2)^(-1)ellが全点で一意に定まり、
Ga4加法と実逆変換からY≅A4×Zがglobal polynomial同型となる。
元の高次項は一次jetの式に寄与しない。β/node/nil等の分母開は除かない。
Zは既約、dimZ=d−9である。

原v02 build_symplectic_sliceの55消去は正重み順で新pivot係数一、他は既消去変数だけである。
該当sourceの全assertionを読解し、global polynomial三角消去であることを確認した。
残134式の保存と再構成の既exact auditにより、Zは実S61全sectionの閉部分として同定される。
旧Q51/Q52 UNIVERSAL_SECTION各§1は全β上の閉A16 graphが各actual normalized成分閉包に等しいと証明済み。
これは単なる全134式への代入からの所属推論ではない。
各C39/C40はd25なので全Zは既約16次元、含まれる閉A16 graphと同次元ゆえ全Z=graph。
よって任意のrank四評価点/実Darboux/実perspective正規化の一つでのgraph testは、
全global T39/T40への所属の必要十分判定である。

独立finite job `jobs/plane_hyperplane_quartic_e6_global_parent_review.json` は実行前gate PASS。
`classification_quartic_second_graph_review.py --e6-global-parent` はproducer/S61算術/CASをimportせず、
E6全六磁気sourceの同一canonical点、rank四点Pf=−75/376、実shift/Darboux/perspective、
全matrixと逆、全150の両向きEuler quotient、全134式と全350Jacobiを再構成した。
全45residualは両graphそれぞれ15本非零で既診断と完全一致。0.367秒でPASS。
旧十FunctionDef ASTは不変、旧全sourcebytes及び実行jobは別snapshotに保存した。
初回の空白欠落Python警告は数学ASTを変えず空白だけ訂正し、旧実行bytesはreport内に保持した。

したがってE6のcanonical実軌道は全T39とT40の外にある。
全source閉包の零magnetic invariant境界まで外と主張しない。
未知親の存在構成、全成分次元や新しい計数はこの補題からは決めない。
A5又は他源への新graph計算もこのreceiptには含めない。
