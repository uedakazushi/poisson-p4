# 918 source の factored rightsolve：独立監査

主報告 `8741d10b...` に対する独立報告 `33d2adb5...` を固定した。
literal A、全 G2 の係数、ordinary quadratic 44順序と ell² の除外、
入力44*target+pair、出力9*source+variable を独立に確認した。
162個の row/column singleton 単元を literal 行列から再確認し、
forward RHS と reverse source の DAG node を一つずつ照合した。
残る28角 block は既存独立証明の非零 determinant と同じ固定行・列である。
その adjugate は展開しない。

保存した918個の source witness は、単元 block の selected original coordinate を与える。
それを「解いた結果」として信用するのでなく、各場合に全28方程式を直接検査した。
2349 scalar node は線形で topological、block の RHS は block 自身に依存しない。
そのため固定単元 block の一意性から、witness は真に同じ factored map R の値となる。
全918出力、全918の1144座標 F R F=F、全252の A R_A A=A、保存 fingerprint が一致した。

符号付き検算の前に、固定 source fixture について正上界134bitを独立に得た。
演算前403bit上界内で逆順の和を使う検算を行い、実 reduced 部分値は最大40bit。
主計算の44bitとの差は加算順序によるもので、任意 RHS に対する上界ではない。
私の検算側では、正 support 上界1711と実 support1253の同一視、slot helper の return、
主計算の44bitとの厳密一致という三つの誤った条件を保存して修正した。
最後の条件は全数学恒等式が通った後の metric 検査だけの失敗である。
主証拠は変更していない。

独立実装の全所有量上界は57955、主実装の66378とは別の allocation である。
report 内 RAW_read counter は継承 hash path を除いた部分計数なので、別の
RAW-accounting-scope-companion が全 unique 入力・出力・失敗記録を含む正確な上界を固定する。
いずれも150000所有量・32MiB・60秒・2GB内で、dense inverse は作っていない。

F R F=F から R は真の像に属する RHS を持ち上げる。
alpha が既知242次元 ker F を消すなら、v-RFv∈ker F より alpha R F=alpha。
これは形式的な演算子の結論であり、実 alpha、Kfull、rho、b1、Z2、
高次応答・support・任意数値 RHS の高さを評価又は保証したものではない。
再現は `code/ambient_B18_ell2_T2_full918_rightsolve_source_review_20261105_v63.py`。
