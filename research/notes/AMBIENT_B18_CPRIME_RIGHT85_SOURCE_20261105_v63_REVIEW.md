# Cprime の一次 syzygy85：独立 SOURCE・保存量監査

2026-11-05。実際の V1 の固定24本を選び、その九一次倍216列を全135二次座標で独立に再構成する。既存 source 証明書の131 normalized pivot と各 source 組合せを掛け戻し、異なる leading row と選択元列への三角変換から rank>=131 を確認する。85 kernel の全135成分消滅と free-tail I85から rank<=131、したがって全一次 syzygy kernel がこの85本で尽くされる。

各 syzygy は文字通り元 V の24本への homogeneous 一次 source 係数であり、V0=0、V1*u1=0なので `V*u1` の三次 initial は `V2*u1` である。高次に変更された pivot 列でも、constant kernel を補正した列でもない。V3・新しい gauge・physical source3を必要としない。

既知三次商の残余 dual は pure r,ell cubic 五項だけに支えを持つ。従って source 係数に z があれば寄与は零であり、V2の nonplane 成分も読まない。この事実を使って85本の folded adjoint を独立に再構成し、463係数・40 active row・115実 cached leaf・最大25項を一致させる。値との内積はまだ計算しない。

adjoint23bit、必要 cached leaf14bitを直接確認し、最大25項の有理数和の930bit上界を得る。固定 source/kernel と入力、85本の最大12項 plane representative、出力を含む22000 numeric envelopeは60000 hard cap以内である。worker の import は独立 `--help` だけで確認し、main を実行しない。60秒・2GB・32MiB制限を維持する。

非零値が一つあれば真の三次 initial が既知 rank26を27へ上げ、既存の上界9・下界9により長い指数e=2となる。全85値が零でも、残る19本の constant-kernel 三次 source はこの試験の外なので、partial negative にとどまる。

監査側の初回 full-hash 転記誤りと、保存 reserve を過大に足した22205という未達 counter は failure-history に保存した。後者は60000 hard capを越えた実計算ではない。source evidence7920、旧二次入力760、metadataと396 leaf・85×25出力・2000 scratch・6000 serialization/controllerを含めた保守上界20205が22000以内であることを再確認した。著者の source・cap・数式は変更しない。
