# 固定12一次切断：未完のsolverと保存済み構造

[証明・正確な限界](../../research/notes/Q_THROUGH3_LINEAR_SECTION_20261103_v63.md)。

整数seed20261103、pivot行列式1の12一次式で35変数から23変数へ写した。全74式のGBはsignal9で最初のcheckpoint前に終了し、全切断の零次元性・次元上界は取得していない。原因は特定しておらず、観測RSSが小さいことから3GiBメモリ不足とは断定しない。

次の軽量な有限証明だけを採用する。

- 既知13平面との交わりはV=Z_(32003)上で正確にV[t]/t³。最大minor46494とcubic値modp11541は単元。
- 全74投影式と有限体行列から、固定切断のHilbert関数は次数3まで(1,23,219,1073)。QQのHilbert等号は主張しない。
- 101個の線形二次syzygyから、特殊ファイバーの四次Hilbert関数は102以上。四次Macaulay全rank証明を排除するが、低次数GBによるArtinian証明は排除しない。
- dense五次Macaulay行列は87.7GBであり実行していない。既存のR35上界22は変更していない。

重いsolverを繰り返さずに軽量証明だけを再現できる。

```sh
source /workspace/cloud-setup/activate.sh
python code/q_through3_section_known_plane_20261103_v63.py --section reports/Q_THROUGH3_LINEAR_SECTION_20261103_v63 --output work/q-section-plane-fresh
python code/q_through3_section_low_ranks_20261103_v63.py --section reports/Q_THROUGH3_LINEAR_SECTION_20261103_v63 --output work/q-section-ranks-fresh
python code/q_through3_section_macaulay_assessment_20261103_v63.py --ranks work/q-section-ranks-fresh --output work/q-section-cost-fresh
```

既知平面0.006秒、tensor投影・rank計3.28秒。全74有限体係数は `low-ranks/projected-polynomials-and-ranks.txt.gz` に保存し、解凍後hashもmanifestにある。全線形形式・solverソース・途中ログ・メモリ観測も保存している。
