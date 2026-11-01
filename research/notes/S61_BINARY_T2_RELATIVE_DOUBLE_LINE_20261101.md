# T2接続の全九変数相対二重線

2026年11月1日。一次項を

    H9={u4=u5=u12=u13=0,5u14−3u49=0,
        5u35+u50=0,u51−5u26=0}, u49≠0

に制限し、Q0=xw∂y∧∂zとする。旧u30=1 chartの残り44 q座標を0から番号付けする。
G=z²wに対応する方向をδとし、非零成分

    n0=5u49, n19=8u26, n24=16u35/3,
    n27=16u49/5, n34=u49

を持つ方向nを取る。このときq=q0+sδ+εn,t=0,ε²=0は全九個の自由一次変数について
全135式を恒等的に満たす。δ32=−2、δ34=0、n32=0、n34=u49なので、
s,εはchart座標から復元できる。
従ってH9上の閉二重線で、その完成環はA[[s,ε]]/(ε²)、A上flatである。

既証明のT2指定点のellはz−15w−23xである。
Darboux shear

    x↦x, y↦y−23w, z↦z+23x+15w, w↦w

はQ0を保存し、ellをzへ移す。
primitive一次項の有限pullbackを全成分で行い、この点がH9へ移ることを確認する。
七つの法線方向へのSp4微分のminorが非零なので、Sp4×H9→W16はその点でsmoothである。
全normalized一次座標・minor・行列をreportに保存する。

相対全135式の完成環S_relからこの二重線D_relへの全射がある。
D_relがβ底A上flatなので、その核N_relの特殊化は固定βの実在二重線への核N0に一致する。
ただしT2点のN0はtorsion-freeと証明されていない。
固定点の初期idealには三次の有限長部分があり、実際のR3-torsionは有限長高々一と分かるだけである。
この補題だけで相対t²=0やambient指定点での冪零性を結論しない。
相対torsionの有限性と飽和境界の次元を別に制御する必要がある。

    .venv/bin/python code/s61_binary_t2_relative_double_line_20261101.py --output work/t2-double-line-recheck.json
