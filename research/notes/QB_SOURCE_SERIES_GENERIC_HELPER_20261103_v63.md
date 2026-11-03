# Q_bのgeneric QQ(b)三source jet helper（author補助）

2026-11-03。凍結済みpoint-QQ helperを変更せず、新しいQQ(b)版を保存した。
これはgeneric cubic gateの**author contribution**であり、そのgateの独立監査ではない。

## APIと係数環

    from qb_source_series_generic_helper_20261103_v63 import make_qb_source

    phi = make_qb_source('CB1', order=3)
    K = phi.domain         #SymPy QQ.frac_field(b)
    b = phi.symbol         #SymPy Symbol
    bg = phi.generator     #FracElement
    jets = phi.zero_jets() #29 offset arrays, degree0..3
    jets[0][1] = bg
    output = phi.evaluate(jets)

kindはCB1、CB2、CA。orderは0..3で、2と3を明示的に検査した。
`parameter=b`も受け付けるが、このfactoryはgeneric QQ(b)であり数値specialization用ではない。
入力は**通常のt冪係数のoffset**で、全定数項は零。
出力は全150物理係数のdegree0..order配列で、全て同じQQ(b)のFracElementである。
入力はFracElementを再parseせずに受け取る。
`.jacobian()`、`.point150`、`.parameter_order`、`.flag_positions`はpoint helperと同じ規約である。

CBの順序はf3,H16,mu,flag9、CAはfiber20,flag9。
全source基点、G_i(b)、H_i(b)は凍結済み実多項式self-pair certificateから取り、有限標本で補間しない。
CBではg(t)=mu(t)*f(t)の全mixed項を保持する。
旗はexp(sum E_ij*s_ij(t))であり、各s_ij(t)の非線形係数も保持する。
opposite flag patternのA^4=0によりgroup exponentialは正確に有限である。

## 二つの正確な評価法

default `.evaluate()` と `.evaluate_physical()` は、5×5外積と逆二次monomialの係数convolutionを行う。
`.evaluate_lie()` は同じ実群作用をexp(L_sum)の係数で評価する別法である。
三source全てで、dense generic QQ(b)三次jetの全150×4出力が両法で厳密に一致した。
二次評価は三次評価の正確なtruncationと一致する。
b=1へspecializeした値は、以前のpoint helperの全150×4出力とも一致する。
specializationは補助checkであり、generic equalityの証明はQQ(b)そのものの両評価法の一致である。

helperの基点・群逆行列が持つ分母は、非零QQ定数とbの冪だけ。
追加のpoleは呼出し側が渡したjet係数の分母からしか生じない。
特定のcubic diagramでb²−3等をinvertするかどうかは呼出し側の証明scopeであり、helperが隠れて追加しない。

## 実行時間と限定

dense generic三次一回はphysical約1.0〜1.3秒、Lie約1.6〜2.3秒だったため、physicalをdefaultにした。
point-QQ版より重く、144回のdense generic評価は120秒の上限に近づく。
呼出し側は許可された一回の有界runと途中保存を維持し、失敗時に無制限のsymbolic展開へ戻らない。

このhelperは87/58拘束frameや逆行列を計算せず、cubic defect、121syzygy、全D局所型の結果を主張しない。
小さなsource評価器の実装とexact checksだけを保存する。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/qb_source_series_generic_helper_20261103_v63.py --output <new-output-directory>

reportと全generic check jets:
`reports/QB_SOURCE_SERIES_GENERIC_HELPER_20261103_v63/`。
