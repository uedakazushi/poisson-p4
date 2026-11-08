# Regular nil HQQ 全 T20 衝突：同じ全変数曲線の有限 readback

2026-11-08、PRIMARY E。root の全六式純手確認と明示 GO の後、
凍結 note `TWO_QUADRICS_REGULAR_NIL_WHOLE_T20_CONFLUENCE.md` と
job `two_quadrics_regular_nil_t20_confluence.json` だけを一回検証した。
fresh gate は input binding PASS。数学検算はその後の一回で **ALL PASS**、
0.053948 秒、最大 RSS13,184KiB、予算45秒／768MiB。

標準ライブラリ Fraction と既存独立 Laurent helper の四則演算だけを使った。
producer／CAS／SymPy／native 又は shared classification runner は import していない。
全 \(e,d,\sigma\) を残し、代入点や新しい operator、rank、minor、profile、jet は作っていない。

確認したのは次の同一有限恒等式である。

* true completed T20 JSON の全 correction を読み、\(t_0=p_1,t_1=-p_1\) による period pole の exact 相殺。
* 元 raw source の全十括弧から fixed GL と h-Casimir Euler gauge を再構成し、全 \(eK+\sigma e^2S\) と一致。
* \(\mathcal L_NK=-2P_0\)、\(\mathcal L_NS=J/2\)、二回目の Lie 零と actual pushforward の符号。
* nil flow、二 translation、fixed frame の全行列と全逆の左右積、total determinant \(1/12\)。
* 表示した六 \(Q_e\) と actual Laurent pushforward の全百五十空間係数、その全 parameter 係数の一致。
  最終式に負 e 冪は無く、全 e degree は二以下。
* 全 raw Euler 補正の正逆、全百二十六 signed Euler-Pfaffian の変換則、全十 strict Jacobi 恒等零。
* 全中心 \(Q_0=P(3d^2/2,\sigma d/4)\)、その cofactor N、actual quintic の完全一致。

実行 inline program 全文、全五行五列 target tensor／total frame／逆、Euler 場、全 actual quintic、
全 input SHA は
`reports/CLASSIFICATION_TWO_QUADRICS_REGULAR_NIL_T20_CONFLUENCE_20261108.json`
（SHA256 `1db8172013c8aff766a00f33877be588d0aaca5578c251b483b1b8ddcb405834`）に保存した。
code を新しい恒久 wrapper にせず、実行全文そのものを immutable report に埋め込んだ。
helper の transitive `classification_quadratic_cofactor_review.py` SHA も明示している。

全 A 非零、全 B と C／旧 \(N^2\) 係数への戻しは、凍結純 note の actual 正逆と閉性を用いる。
この一回は、その命題の別の点を探したものではない。
root の独立純監査・正式採録 receipt と、今回の同作者による全係数 readback は区別する。
新 component 数は零、他の未認識 nil 軌道を閉じたとはしない。
旧 proof／job／manifest と shared 0b7f／7619 code は全て不変である。
