# 既存低不変式の既知分母Hを正確に約分する計算

2026-11-04。受理済みのF4、p0,...,p4、N12の整数係数だけを使い、既知分母
H=147456Y−(3072W−1)²の因子をQQ多項式除算で取り除く。一回60秒・2GB。
各出力で元の整数多項式への再乗算を検証する。新点・方程式探索・次数上昇は使わない。

実際の約分数と項数は同名reports directoryのreport.jsonへ記録する。
これは全不変式環やNEW方程式、target函数体indexの決定を主張しない。
出力の係数は分子・正分母をそれぞれ十進文字列で保存する。

再現は新しいoutput先を指定して
`.venv/bin/python code/q_new_family_low_invariant_pole_cancellation_20261104_v63.py --output work/low-pole-replay`
とする。既存七写像の全係数を変更せず、商とH冪の積がそれらに等しいことを検証する。
