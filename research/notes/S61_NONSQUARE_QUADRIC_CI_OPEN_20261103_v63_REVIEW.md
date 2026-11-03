# 明示的 CI-open の独立監査

2026-11-03。`S61_NONSQUARE_QUADRIC_CI_OPEN_20261103_v63` の actual source/resultant と幾何学的 height 論証を監査する。**R_vertex は追加の未知多項式条件であり、全 core の既知単元ではない。**

元804項、四座標 frame、C1/C2、全三つの determinantal minor と binary D は、既に独立受理済みの binary reduction review に固定される。同じ source の actual pair chart と全自然係数は先行 independent finite review に固定される。全入力 hash と今回の frozen archive を照合する。著者の producer は再実行せず、既存の独立 source 認識を依存として用いる。

pair chart は A=X=W=0、B=−10aPJ*z であり、aPJ が単元なので vertex A=X=W=B=0 は z=0 と同じである。全12 pair 係数から z を含まない三係数を取り、二つの homogeneous binary quadratic の Sylvester4角を作る。独自 Fraction 辞書の24項 determinant は厳密に

`−729*t*R_vertex(t,v,H)/4096000000000000000`

となる。保存された R_vertex の全92係数・総次数9を確認する。binary の先頭係数の退化も homogeneous resultant が扱い、追加の先頭係数の逆元は用いない。非零は vertex P1 上で共通零点がないことと同値である。点(1,23,9)の値はこの chart の非空性だけを確認する。

## 退化した curve を含む幾何学的証明

以下は各幾何学的 characteristic-zero 点上で行う。frame は aPJ のみを分母とし、A,X,W,B が独立座標である。Q1,Q2prime が共通 quadratic factor を持てば比例するが、W=0 上の AB,40XB は比例しない。共通 linear factor は W ではなく、modW で B に比例するので B+kW となる。`W*D=40X*Q1−A*Q2prime` からこの factor は D も割る。しかし B=W=0 上の D は X²係数 `1280a³J²/243` を持つ非零 binary なので矛盾する。よって Q1,Q2prime は互いに素である。

三つの2-minor ideal はこの互いに素な二式を含むので height≥2。一方2×3 maximal-minor の height bound は2であり、height=2となる。Hilbert–Burch により四変数 quotient はCM、projective support C は純一次元である。既約性・reducedness・smoothnessは不要である。

集合として `V(Q1,Q2prime)=V(B,W)∪V(全三minor)`。W≠0 では上の恒等式で D=0、W=0 では B=0 または A=X=0 であることから従う。後者は六変数空間で C を base、二変数 S,U の P1 を vertex とする cone である。

R_vertex≠0 とすると、四式の共通零点 Y⊂P5 は vertex P1 を避ける。P3=V(B,W) 内の閉部分 Y∩P3 が次元2以上なら、projective dimension theorem によりこの P1 と交わるので、次元≤1である。

cone 上の Y の部分も projective な閉集合で vertex を避ける。vertex からの射影は C への morphism となり、source が projective、target が separated なので proper である。各幾何 fiber は ruling P2 の閉部分であり、その vertex line P1 と交わらない。正次元なら projective dimension theorem に反するので全 fiber は零次元以下である。従って proper quasi-finite、特に source の次元は C の次元1以下となる。この論証に非被約 base/fiber の reducedness は不要で、support の dimension を使う。

両部分から dimY≤1、六変数の四生成元 ideal は height4となる。多項式環はCMなので四式は regular sequence4である。別の gcd-Σ 排除を追加する必要はない。

## 元の二次式への帰結と限界

FULLU と D(R_vertex) の共通部分では ell chart が可逆なので、raw CI は元の四式 mod ell のCIと同値である。ZW product19 の独立レビューを合わせて colon exact sequence が各点で成立し、HF4=48+16=64、元の light quartic products の rank は146となる。

同じ列を次数2に使うと、六変数CI4の HF2=21−4=17、ZW五変数商の HF1=5 なので HF2(S/I)=22である。従って元六二次式はこの開の全点で独立（28−22=6）。heavy×六二次式の列は light 四次と disjoint なので weighted quadratic-products の rank は152となる。

R_vertex=0 上のCI、全点下界146、full cubic57、残る quartic quotient の全射・長65・ambient除去は未解決である。全 core の既知単元として R_vertex を扱うこと、標本点だけから上の全点結論を出すことはしない。

再現（新しい出力先）:

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_quadric_CI_open_review_20261103_v63.py \
  --output work/marked-square-v63/CI-open-independent-new
```

60秒・2GB以内で、GB、新しい体、有限体、巨大四次行列を起動しない。
