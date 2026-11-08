# simple-factor spectral input の polebudget 文章の訂正

2026-11-08。入力 `inputs/nil31_short_simple_h_spectral_curve_20261108.json`
の `polebudget_proof[0]` だけに、旧 frame の `diag(1,e^2,e,1,e)` という文章が残っている。
新 frame では正しく **diag(1,e^2,e,e,1)** である。
同 input の実データ `one_parameter_weights=[0,2,1,1,0]`、constant M/inverse、
spectralframe、executed program と全 operator/endpoint matrix は正しい新値を使っている。
新旧とも quadratic h-Casimir weight range[-3,3]、Laurent bound[-4,2]なので、
保存した missing-monomial obstruction と全150 target transfer の数学は変わらない。
executed input/raw report/成功した負の証明の bytes は訂正せず保持する。
この誤記を old-H frame で演算した証拠、又は新 endpoint matrix の source として使用しない。
