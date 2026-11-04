# h22 の三純行・固定三列の実小行列式

受理された source/height preflight の候補を一つだけ実際に計算する。源は実135ラベル18,39,68、標的は元32列11,13,25。成分は h22 P_i[j]-P_i[22]Q9[j]。h22 は既存十四項係数であり、他のQ9尾を零とはしない。三源は literal Hzero なので、Hの随伴行列・B分母・未知の逆因子を導入しない。

既存108 scalar node の入力を source binding とともに抽出し、a の Laurent 単項式と QQ 係数を扱う。未知パラメータでの除算や一般 gcd は使わない。全九成分、三つの2×2 cofactor、cofactor 展開と六置換展開が与える同じ実3行列式を記録する。全108中間値も保存する。実値・零非零の結論は report.json に固定する。

資源ガードは入力1120項、全scalar値の support 上界、九成分4630項、行列式10080項、cofactor・六置換の作業領域を合計する。native は全生存多項式の実項数と入力reserveを毎演算で数え、150000を超えない。実行全体60秒、2GB、成果物32MiBを守る。

零ならそのまま保存し、別の列を探索しない。非零でもこの三純行ブロックの一般的独立性を与えるだけであり、六連成行を含む全九核の独立性、全 FULLU の三次57、四次239は証明しない。既存Delta5/h22チャートの補集合も残す。高次因数分解や単元剥離はここで自動実行しない。

再現:

    source /workspace/cloud-setup/activate.sh
    .venv/bin/python code/s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.py \
      --output reports/S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63

全入力SHA、抽出した108根と原子、native完全入力、全中間値、九成分、二つの行列式表示を同名報告に保存する。
