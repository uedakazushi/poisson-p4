# θ/H の表現を保持した全24次基底への修正：事前検査

2026-11-04。C/D trace6だけを使った一つの形式平面検査はrank8で停止した。旧Inv₁₂基底が使うのは、4×4の二つのadjoint C,Dではなく、

\[
 V=(\Omega\rho(A_j)\beta)_j,\quad
 \mu_j=\tfrac12\beta^T\Omega\rho(A_j)\beta,\quad
 \theta=\rho(\operatorname{Gram}^{-1}\mu),\quad
 H=\theta+\Omega^{-1}V\operatorname{Gram}^{-1}V^T
\]

という**16×16の作用行列とHessian**である。Hがsp₄の4次元表現に属するとは仮定しない。旧`DEGREE12_14_INVARIANT_EXCLUSION`の全9次元基底はこのθ/Hのtrace6である。

同じ三つのambient vectorと同じ形式平面を保持し、保存済みC/D rank8空間へ旧14 trace6を加える案を検討した。必要なら一つの旧wordだけで9次元を補い、F₄Inv₁₂+Sym²Inv₈のproduct rank19を確認する。16次の追加成分には保存済み36 C/D trace8を先に試し、不足する場合に限り同じ既存のθ/H trace8を使う。最後に70列から67列の独立なものを選ぶ。別平面・追加のactual NEW25点は使わない。

形式平面のθ/H各成分は二次、係数数は高々6。全trace6は長さ2の4個、長さ3の8個の行列積と14個のtrace pairingで計算でき、scalar polynomial productは52,736回以下。全trace8を追加する場合はさらに74,752回以下。30個のprefix matrixを全て保持しても、degree8までの係数上界は345,600個のmod32003 residueであり、2GB上限には十分余裕がある。実時間の保証に置き換えず、native一回60秒・2GB・32MiBで停止する案とした。

将来のsource展開については、既存のβ整数numeratorと定数行列だけから正のℓ¹係数majorantを計算した。**新しい不変式の係数展開やrankはこのpreflightでは計算していない。** θとHの全quadratic係数に対し共通整数scale20が十分。全wordでの上界は次の通り。

|source word|β次数|B取消し後の項数上界|整数係数bit上界|compact JSON bytes上界|
|---|---:|---:|---:|---:|
|旧θ/H trace6|12|27418|506|5291674|
|旧θ/H trace8|16|80125|674|19470375|

B分母次数d/2以下という旧preflightのweight論証は表現行列の大きさに依存せず、これらにも適用できる。H分母はd。source上で全prefixを同時展開せず、選択wordを一つずつhalf-word二枚と一時領域だけで処理する必要がある。runtimeやmemoryの完走保証は未取得。

CD trace6二本・旧trace6一本・CD trace8一本の最悪raw保存上界33,620,451 bytesは、32MiB=33,554,432を僅かに超える。最後も旧trace8なら34,742,201 bytes。従って後続source展開には、実際の選択wordについての上界再確認と監視付きcap、または保存表現の検討が必要で、成功を事前に保証しない。

保存は`code/q_new_family_degree24_theta_H_repair_preflight_20261104_v63.py`と
`reports/Q_NEW_FAMILY_DEGREE24_THETA_H_REPAIR_PREFLIGHT_20261104_v63/report.json`。全旧word、定数入力hash、representationの区別、演算数と係数majorantを保持する。bound metadata計算は約2.82秒。

この後に別の明示的budgetで実行した同一平面の修正版は、`Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63`に保存する。本preflight自体はbasis・source関係式を主張しない。
