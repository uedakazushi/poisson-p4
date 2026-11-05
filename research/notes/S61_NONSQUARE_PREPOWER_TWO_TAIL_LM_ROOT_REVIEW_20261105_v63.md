# 二tailの六因子base環における非割切り：独立確認

元4500項N0と、全整数掛け戻しを済ませた2317項N2のleading monomialを別stream parserで求め、lex a>e>dでそれぞれ(16,10,18)、(12,11,14)と確認した。e指数の差は−1なのでN2はN0を普通多項式環で割らない。

完全なbase環はQQ[a,e,d,1/(a*d*F*B*C*J)]である。六因子の順序を逆にして積を再構成し、元FULLUの全49係数と一致した。a,dは変数、F,Cはaの一次式で定数leading係数500,100、Jはeの一次式で定数leading係数−6。Bはdの一次式で、定数項60e−40の唯一の因子3e−2は先頭係数を割らない（e=2/3の剰余は50a）。従って全六因子は既約で、literal scalar normalizationによって互いに非同伴である。

N2はa,dを割り切る最小指数が零で、F,Jは既存のv(Q1)=0、Cはexact vC(T1)=4、Bはexact vB(Q1)=4と全B四乗恒等式から、六因子全てと互いに素である。N2がs*N0を割りsがこの六因子の積なら、UFDでsを除けるためN2|N0となって矛盾する。従ってこのbase環で(Q2)は(Q0,Q2)に真に含まれる。

これは二つのtailのscheme idealだけの結論である。D12、Tや未知の別因子をさらに反転した環、strict open、radicalの相違、実bad Fittingや空補集合は主張しない。固定有限体gcdが正であったことから標数零の新因子も主張しない。

独立報告 `reports/S61_NONSQUARE_PREPOWER_TWO_TAIL_LM_ROOT_REVIEW_20261105_v63/report.json` は全入力SHA、逆順49係数積、stream LM、受理済valuation証拠とsource codeを保存する。
