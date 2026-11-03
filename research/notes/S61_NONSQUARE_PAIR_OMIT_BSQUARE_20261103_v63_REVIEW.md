# B²を除いた pair の独立算術監査

2026-11-03。`S61_NONSQUARE_PAIR_OMIT_BSQUARE_20261103_v63` のうち、**全 aPJ-open 上で actual2×5係数行列が rank2** という算術定理のみを受理する。著者 note §2–3 の syzygy/rank78 への適用は受理しない。独立した gap note を併記する。

元804項、既知 unit 正規化、全自然係数、f=0以外の C/Q23 候補の証明は、先行する独立 whole-pair review と同じ frozen source である。全入力 hash と今回の archive を照合する。候補決定に用いた列2,3,4とu²列5の三つの minor は、B²列0を削除してもすべて残る。C枝の除外、Q23内の Hbar,tbar の決定には列0を用いない。

Q23を法として、actual 第1列とu²列の minor を全係数から作り、Hbar,tbarで独自 Fraction 演算により評価する。保存された inverse との積が1であることを一変数長除法で確認する。Q23の既約性・reducedness・体構造は仮定しない。

f=0 の二値 v=−3,−4/3 では、Hを指定する前に `row0_col3−(4/15)row0_col4` がそれぞれ `(3/625)t`、`(41/20000)t` となることを確認する。tは既知unitなので第1行は残した列だけでも非零。第1行u²係数は零であり、rank failureなら第2行u²係数も零となる。その係数はHに関して非零一次なので、H=9,7/3 がそれぞれ強制される。

両枝で残る列1〜4の全六 minor を actual 係数から再構成する。t,Pと非零QQ定数だけの保存因子を確認し、六つの primitive minor に対する literal QQ Bezout sum=1を独自 Fraction 辞書で確認する。従って両枝も排除される。

以上で全幾何学的点でのrank failureが空となり、aPJ局所化の有限型係数環上で maximal-minor ideal はunit、cokernelは零となる。これはB²を除いた stronger pair criterionであり、ternary pairが共通lineを持たないこととは異なる。

**未受理部分。** 著者 note §2 の B*S/B*U 係数比較は scroll側linear cofactorの寄与を落としている。実際に S*Q1 は A*B*S/6075 を含む。従って2×5 rank2から五式のlinear-syzygy空間が二次元とする推論は未証明であり、§3のG1/G2からraw rank78への適用も受理しない。著者は指摘を認め、この適用を撤回した。元のfrozen bytesを変更せず、この監査とgap noteで範囲を訂正する。

既存の全域product41は別の**五式の定数関係**の排除を用いており、今回の五式のlinear-syzygy推論を使わない。先行のQQ41・CI-openの受理範囲は変わらない。

再現（新しい出力先）:

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_pair_omit_Bsquare_review_20261103_v63.py \
  --output work/marked-square-v63/omit-Bsquare-independent-new
```

公開reportの `all_passed_within_explicit_arithmetic_scope` は上記算術だけを表す。二つのsyzygy/rank78の accepted flag は明示的に false とする。新しい探索、GB、体、prime、重い source 再展開は行わない。
