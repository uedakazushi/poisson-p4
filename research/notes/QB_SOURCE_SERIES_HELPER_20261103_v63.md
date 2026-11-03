# Q_bのpoint-QQ三source series evaluator（author補助code）

2026-11-03。CB1、CB2、CAの実29座標sourceを、Q_bの指定非零QQ点で評価する小さな補助code。
高次probeの**author contribution**であり、そのprobeの独立監査ではない。
凍結済み研究packetを編集せず、新しいhelperとして保存する。

## API

    from fractions import Fraction
    from qb_source_series_helper_20261103_v63 import make_qb_source

    phi = make_qb_source('CB1', parameter=Fraction(1), order=3)
    jets = phi.zero_jets()       #29 rows, degree0..3
    jets[0][1] = Fraction(1)     #source offset coefficient of t
    output = phi.evaluate(jets)  #150 rows, degree0..3, exact Fraction

kindはCB1、CB2、CA。orderは0、1、2、3。
入力は凍結済み実source基点からの**offset**で、全29行の定数項は零でなければならない。
入力・出力とも通常のt冪の係数であり、k階微分やk!倍ではない。
短い行は零でpadする。指定orderを超える行と非零定数項は拒否する。

- CBの順序はbinary f3、Hamiltonian16、ratio mu、全flag9。
- CAの順序は元のfixed-fiber20、全flag9。
- `.parameter_order`、`.flag_positions`、`.point150`で正確な対応を取得できる。
- `.jacobian()`は同じ順序の物理150×29 QQ行列を返す。
- `.evaluate_lie(jets)`は同じ実群exponentialの別の係数評価器であり、cross-checkや計算の比較用。

## 実sourceの保持

CBは3binaryと16cubic係数にmuを加えたsourceであり、g(t)=mu(t)*f(t)を毎回convolutionする。
全mixed次数2、3を含み、自由な三g係数でratio一座標を置き換えない。
CAは正確なlinear20 sourceを使う。
群作用は

    exp(sum_{flag9} E_ij*s_ij(t))

であり、ordered product of exponentialsではない。
flag patternには長さ4の非零matrix pathが存在しないためA(t)^4=0。
従ってI+A+A²/2+A³/6は実際の有限matrix exponentialであり、非線形なs_ij(t)も保持される。
標準CB frameから元の物理Q_bへ戻す実G_i^-1も含む。

物理評価器は外積5×5行列と逆二次monomial代入を、tの係数配列の有限convolutionで計算する。
全計算は`Fraction`上であり、浮動小数点やmodp評価を使わない。
入力の正確なsourceは親の実Q_b certificateと全CA/CB基底へhashで結び付ける。

## 小さな検証とscope

各三sourceで全29単位一次方向を評価し、全150×29 Jacobianと照合した。
全source係数を動かすdense三次jetについて、物理座標評価器と独立のexp(L_A(t))係数評価器の全150×4出力が一致する。
order2評価がorder3評価の正確なtruncationとなることも確認した。
CBではdelta fとdelta muを同時に動かし、非零のmixed二次係数を明示的に検査した。
零offsetは全150で凍結済みQ_bの点へ戻り、高次係数は零となる。

degree3のdense一評価は約0.04〜0.08秒。
これはsource評価器の検証であり、87拘束のsource inversion、121syzygy、cubic gateの結果は主張しない。
その計算を行うagentのauthor補助として提供する。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/qb_source_series_helper_20261103_v63.py --output <new-output-directory>

結果: `reports/QB_SOURCE_SERIES_HELPER_20261103_v63/report.json`。
