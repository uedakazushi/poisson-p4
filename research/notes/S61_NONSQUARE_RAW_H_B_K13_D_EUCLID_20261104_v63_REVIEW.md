# Ka13 全体の d-s common-root cover の独立監査

凍結 D_EUCLID の四つの LC 単元証明、十五 polynomial nodes、source cofactor SLP と最終
三除算・Bézout 恒等式を独立再生する。作者の商環 arithmetic や Euclid worker を使わず、
固定 table a^0,...,a^24 の reduction を整数恒等式で検証した別実装を用いる。
新たな係数逆元・GCD・solve・GB・点・素数は探索しない。

H4 は独立受理済み実 monic coefficient に結び、g,k は元の actual source minors15,20の
147項を Ka13 で直接評価し、保存された H4 quotient/remainder と全係数比較する。
各保存 LC について primitive/content と整数 S*p+T*g13=R≠0 を検証し、保存 inverse
との積が1になることを確認する。従って体性・既約性・reducedness は不要である。

source/単元正規化/monic division の十五 node を、source cofactor SLP を含めて確認する。
最終 node10 は monic d-s で、三 source inputs はすべてこれで割り切れ、逆に d-s は
三入力の明示 Bézout combination である。よって Ka13[d] 上の ideal は正確に
(H4,gR,kR)=(d-s)。元の g,k への変換は保存された H4 quotient を係数に吸収するだけで、
別の逆元を使わない。

これは必要 common-root cover の scheme 等号である。残る actual B minors を消したとは
主張せず、actual Bbad の存在・空性、full cubic57・quartic239 を結論しない。
結果・時間・入力 hash は report と independent-proof に保存する。
