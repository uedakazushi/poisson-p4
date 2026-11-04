# Companion随伴列からのe復元と13次projection代数

2026-11-04。一個の保存V成分の単元証明と全16個のEV=tVだけを計算する。e-Euclid、新norm、finiteNF Gaussianは使わない。

H=h1のmonic16基底でe乗法companion Eを作る。保存Vnormal mod g13の各成分をnumeric分母で整数化し、一成分のprimitive numeratorがg13と互いに素であることをliteral整数Bézoutで証明する。Vを成分別にscaleし直さず、inverse構築の補助numeric contentだけを記録する。

t=(EV)_j/V_jをg13上で計算し、全16個のEV=tVを完全な整数divisibility恒等式として保存する。未知のa多項式を証拠なしに反転しない。

companionにはh1(E)=0が成り立つ。旧norm/随伴列の受理済み恒等式からNkV=0、mixed W証拠からNm1V=0。全EV=tVよりp(E)V=p(t)Vとなり、証明済みunitV_jを使えばh1(t)=k1(t)=m1(t)=0。従ってT=QQ[a,a^-1,e]/(h1,k1,m1)からKa13=QQ[a]/g13への全射がある。既存Fitting長さ上界13とdimQQ(Ka13)=13から同型となる。この論証はg13のirreducibilityやkernel-rankを仮定しない。

これは明示triple projection代数のexact length13、monic H4 coverのlength52であり、actual Bbadが非空だと示したものではない。source全22条件とFULLUをさらに課す必要がある。G0側は別に保持する。

e復元後はoriginal22 source、FULLU product、monicH4、代入a/eの全定義とinput degree/heightを束縛し、4×88 Ka13 modelの準備costだけを記録する。そのnormal form係数、target membership、Gaussianは次の予算に残す。

## 今回の確定結果

承認された一回はnative0.4166秒、全体2.0467秒で完了した。選択pivotはV15、primitive numerator9277bit。unit Bézout S/Tは13/12項で、全16個のEV=tVと共に整数恒等式を保存・再検査した。復元t(a)は13項、numerator1444bit、denominator1408bitである。従って明示triple projectionはKa13と同型、長さ13、H4-cover長52となる。actual badの実在は結論しない。

全22 original sourceは2963項、degree(a,e,d)=(7,11,7)、FULLUは49項、monicH4は74項でa最小指数-5。将来4×88行列は最大4576個のQQ係数を持つ。しかしtの実heightからconservativeなpseudo-reduction heightを取ると、全matrixのRAW上界は32MiBを超える。この上界を小さいと偽らず、field matrixの自動展開はしない。

代わりにH4係数だけなら52個のQQ slotで、既知g13/tの整数heightから得るRAW boundを保存した。これだけを先に準備すれば、実H4 heightを使って次のmatrix予算を絞れる。今回H4 normal coefficients自体も計算していない。
