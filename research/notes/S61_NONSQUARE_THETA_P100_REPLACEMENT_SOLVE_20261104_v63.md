# P100 の二つの source minors による全点 cubic57 判定

先に凍結した 17×24 行列の既知単元前処理を入力とする。これは P100∩FULLU 全体の
行列であり、既存の generic-rank packet の同じ source rows と columns を保持する。
準備済み A17×17、残る B17×7 に対し、一回の exact fraction-free solve
AX=delta B を実行する。全 119 cofactors を local checkpoint に保存してから、
残る R46 と互いに素な最初の replacement entry だけを探す。

実行時に delta が、保存済み degree447 minor を既知の行・列単元で割った
expected degree134 determinant の非零 QQ 定数倍であることを検証する。さらに
全 408 source-unit identities と選んだ 17-vector の AX=delta b を、native solve
とは別の整数係数積和でも確認する。R46 との互いに素性は整数係数 S,T と非零整数
C を用いる literal identity S*Xij+T*R46=C で証明する。容量上限に入らなければ
cofactor checkpoint を残し、この Bézout 保存を未完として扱う。

全係数恒等式が保存・確認された場合、R46≠0 の点では元の selected minor が単元、
R46=0 の点では選んだ replacement minor が単元となる。Cramer identity により
replacement entry は実際の source minor の非零 QQ 定数倍であり、前処理の
全係数は chart 単元である。従って二つの minors が全 P100∩FULLU を覆う。
既存の source projection theorem により remaining rank17、augmented rank25、
実際の full cubic57 が全点で従う。

実行結果と各検証の成否は report.json と independent-verification.json に記録
する。この判定は full quartic239 や P100 の chart 外、P80 の新規結果を主張しない。
既存の P80 algorithm を同じ上限 60秒・2GB・32MiB で再利用し、K46 Gaussian、
数値点、有限体、GB、別の selected-minor search は行わない。
