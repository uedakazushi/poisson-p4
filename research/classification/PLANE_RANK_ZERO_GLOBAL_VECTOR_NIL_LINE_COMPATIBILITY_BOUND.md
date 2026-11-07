# normalss regular-nil二parameter源の最後のK零線を全supportで閉じる

2026-11-07、主義務G。producer exact characteristic-zero証明。
`PLANE_RANK_ZERO_GLOBAL_VECTOR_TWO_PARAMETER_BOUND.md`のnormalss NILREGについて、
最後に残ったwhole K零線を処理する。

全first symbolはnormalssで`T=tJ+J²`, `V=J`, `K=0`、t∈C。
全normalizationのinverseは既存semisimple first-jet noteにある。
fixed-plane affine symbol baseは十、tを固定した層では九。
full thirty/nine higher source、actual有限shear、full inhomogeneous weight二solutionを保持する。

## 1. generic lineの全binary compatibility

full D1 kernel八、actual shear六なので、全first weightをactual有限shearで二変数u0,u1の
transverse sourceへ直せる。
任意のfixed Π1に対するΠ2 solution fiberはfull D2 kernel三のtranslateである。
full Q2=(1/2)E∧[Π1,Π1]をD2の全cokernelへ射影したcompatibilityの中に

`4*u0²/(t+1)²`、`12*u1²/((t−1)(3t−1))`

が現れる。reportはそれぞれのfull cokernel row132/134と全source basisを保存する。
二式のbinary Sylvester resultantは

`2304/((t−1)²*(t+1)^4*(3t−1)²)`

で非零。
rank、basis minor、kernel分母及びこのresultantの全exception factorを戻すと、
例外は正確にt=±3,±1,±1/3だけとなる。
その外では全reduced first-weight compatibilityはu0=u1=0なのでdimension零。
したがってfull higher dimensionはshear六＋weight二fiber三＝九以下、
ambient imageは`symbol10+higher9+plane6−scale1=24`以下。

## 2. 全六exceptional points

全六値で元のfull350×30/9/shear matricesからkernelとactual transverse sliceを取り直した。
generic sliceを特異点へ代入する仮定は使わない。
各sourceでfull weight二compatibility quadricsのgcdは1であり、polynomial UFDでheight二以上。
unit ideal又はsource空という主張ではない。

| t | nullity D1 / D2 | actual shear rank | transverse weight一変数数 | higher上界 | ambient上界 |
|---|---:|---:|---:|---:|---:|
| 3,−3 | 9 / 3 | 6 | 3 | 10 | 24 |
| 1,−1 | 9 / 3 | 5 | 4 | 10 | 24 |
| 1/3,−1/3 | 8 / 3 | 6 | 2 | 9 | 23 |

特にt=±1ではactual shear rankが五へ落ち、sliceは四変数である。
前noteの未実行課題に記した「全exceptionでshear六、三変数slice」という見積りをここで訂正する。
full再構成ではその一係数を捨てず、height二を四変数の全idealで検証した。

t=3ではrow64の`−u1*u2`とrow132の`u0²/4`はcoprime。
t=−3ではrow131の`−u0*u2`とrow132の`u0²−2u1*u2/3`はcoprime。
t=1ではrows32,33,132の`−2u0u3`, `u3²/2`, `u0²`のgcdが1。
t=−1ではrows32,51,133の`−u2u3`, `3u3²/4`, `−2u1u2−3u2u3/4`のgcdが1。
t=±1/3では両変数のsquareがcokernel quadricsに入り、全compatibilityはoriginだけ。
すべてliteral full sourceのexact式であり、特定higher点の不成立から次元を推測していない。

## 3. 全二parameter源への帰結と証明書

K非零は既存two-parameter full-wall certificateで閉じていた。
本命題でK零のwhole t-lineと全exceptionを処理したので、normalss NILREG(t,k)全源がG非優越となる。
二parameter四族はこれですべて閉じ、normalss MIXHは二十三以下、MIXNは二十二以下、
normalss NILREGは二十四以下、normalnil b零generic mixedは二十二以下。
別表のspecial mixed rowsを戻してもnormalnil全mixed m二源は二十三以下。

exact immutable certificateは`CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_NIL_LINE_20261107.json`。
full operators/kernel/basis、全Q2/cokernel三百五十式、binary resultant、全exception pointのexact gcdを保持する。
reusable runnerのisolated `--global-nil-line` modeで、旧二十FunctionDef ASTを保存した。
gate PASS、240秒/1536MiBのfinite contract、実測4.863秒/67,168KiB。
nonlinear GB又はjet次数追加は使わない。

normalss generic MIX三parameterとnormalnil b零Cartan三parameterは依然別義務。
scalar co-Higgs、pointwise-dependent trace/common field、高mに復活させるsourceも残す。
全Gの完了ではない。
