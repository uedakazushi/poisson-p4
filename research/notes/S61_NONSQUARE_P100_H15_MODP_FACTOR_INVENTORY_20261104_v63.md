# P100 の有限候補 H15 の既存素数での因数一覧

入力は凍結済み `S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63/native.json` の exact_gcd_monic のみである。これは既存の五つの 9 小行列式のイデアルを定めるモニック十五次多項式であり、実四次欠損とは同定しない。

元の QQ 係数の分母を全て検査して既存素数 32003 で還元し、FLINT の nmod_poly_factor を一度だけ呼ぶ。全因子の積と重複度を再構成し、各因子に対して Rabin の既約性判定を独立に検証可能な形で保存する。具体的には全 Frobenius 鎖、最終 x^(p^d)=x、d の素因子 l ごとの gcd(f,x^(p^(d/l))-x)=1 の多項式 Bezout を記録する。

この因数分解は有限代数全体を有限体上で計算する際の資源見積りにのみ用いる。新しい素数、beta 点、zeta の一点での評価は使わず、標数零での因数分解や実四次階数239を結論しない。モニックなので Z_(32003)[z]/H は有限自由階数十五であり、後続の正の証明は全因子・全重複度を含む必要がある。

再現:

    source /workspace/cloud-setup/activate.sh
    .venv/bin/python code/s61_nonsquare_P100_H15_modp_factor_inventory_20261104_v63.py \
      --output reports/S61_NONSQUARE_P100_H15_MODP_FACTOR_INVENTORY_20261104_v63

上限は一回 60 秒、2GB、成果物32MiB、明示係数150000。実行結果、因子次数、重複度、全証明書と入力 SHA256 は同名報告に固定する。後続の source-only preflight と実係数・階数の仕事は分離する。
