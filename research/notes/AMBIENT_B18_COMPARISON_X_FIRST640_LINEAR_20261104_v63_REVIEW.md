# 実X1の最初640列：独立prefix・source4・dual監査

2026-11-04。作者worker/helperは呼ばず、既受理source係数に対する疎dual代数
`QQ[zA0,...,zA6,r]/m²`を独自の辞書積で実装する。新しい物理jetや陰関数求解はない。

実chain operatorの17849 nodeから全2808 scalar rootの祖先を独立抽出し、1302 nodeと
79個の元source関数に一致することを確認する。既存source105/H20の値から、同じCA-F9座標の
定数・一次係数を抽出する。product idealはm²に含まれるため、H20上のfirst coefficientは
この実regular baseのfirst coefficientである。旧792/784 cacheは元node IDで固定する。

35個のordered CA prefix商は微分を使わず、保存source2のmonomial supportから再構築する。
第j商に寄与するのはz_jを含み、後続zを含まないmonomialだけである。
z_jを一個除けば、z_j²の係数はそのままlinear z_j係数になる。factor2は掛からない。
r,sは全prefixで保持し、最後に実ell=0のs=rで制限する。

八つの旧CB/新CA regular g商については、先行独立source4監査を固定した上で、
各族のnumeratorの全15係数に `numerator4=g3*quotient1` を再検査する。
membershipと商原点値零は既存全次数の証明であり、この四次だけからdivisibilityを推定しない。
CBのplane勾配を保存された実full9 adapterでCA-F9へ運び、全八勾配を独立復元する。
CB raw f63 node28539とplane gB node299650を区別し、raw f63のfirstはown source2とadapterから
別に比較する。原chainの最後のmultiplierをnew gAへ置換しない。

同じregular origin cacheを使うoriginal C0 DAGの疎ADで79 atomとchainの全1302祖先を評価する。
P8は割当て済みの方法に従い、全source4/P640を独立再評価した先行QUOTIENT_FIRSTS_QUARTIC
peerの正確なhashに束縛されたfirst640×81 cacheを使う。単位constant target/source消去の後に
残す8行×640列の勾配を実s=r adapterで制限する。全56×640の定数は先行actual X0 peerに一致する。
出力の全409非零QQ係数を確認する。

追加のfresh P9 original-DAG再評価も試みたが、割当てられた79 chain atom closure外の
元CA_trace_p node11で不足するcacheを検出した。この任意の追加経路はfrontierとして保存し、
その微分を推測せず、新sourceを問い合わせず、既に独立source検証済みの上記P8 cache経路を使う。
Pのfresh再構築を本peerの実施内容に含めない。

旧r/gは証明済みの正則germのままであり、零locusへ代入して再除算しない。
native128列X1、Gaussian、Euler wedge rank、Tの最小生成数・長さ、K torsion、元深さは未計算。
専用checker `code/ambient_B18_comparison_X_first640_linear_review_20261104_v63.py` と同名REVIEWに保存する。
