# B18におけるCB正規化sourceの全分岐scheme：独立監査

2026-11-03。CA/Nの監査と同じ独自の150係数Lie作用を使い、CBを別に再構築する。固定旗の22基底を物理的な式から作り、全16 parabolic作用とratioのcharacterを含む19列の相対安定性を確認した。

ratioとbinary係数の全積を残したsource写像を微分する。中心の固定19係数、parabolic orbit10 minor、全断面の可変fiber20列、旗9列を検証した。fiber blockの行列式−864/25とflag8 blockの行列式8/3は断面全体で単元であり、全最大minor idealを全Schur列へ正確に置換できる。全150成分を保存値と照合し、lex順序の別GBで予測idealとの両包含を確認した。

全source ramification cone完成環は `C[[g19,w,u,v]]/(wu,wv,uv,u³)` となる。CAのsource ramificationと抽象的に同型であるが、実targetのCB完成環やCA/CB交叉Dを同一視しない。primary構造とdepthの小単項式計算はCA/N独立監査を再使用する。projective版は滑らかな因子を一つ除く。

検証コード `code/ambient_B18_CB_critical_review_20261103_v63.py` はproducerの作用関数や実行関数を呼ばない。実odd ideal、実Dのnilradical・Ass・depthは未解決である。
