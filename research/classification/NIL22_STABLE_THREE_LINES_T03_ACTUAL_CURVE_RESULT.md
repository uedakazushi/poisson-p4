# Stable three-lines 全六磁気源の T03 actual confluence：有限恒等式結果

2026-11-08。純証明 `NIL22_STABLE_THREE_LINES_T03_ACTUAL_CURVE.md` を変更せず、
その同一全変数曲線の算術結果を別に記録する。状態は **producer 全恒等式 PASS、別作者 review 待ち**。
新成分数は零。stable_line_conic には適用しない。

## 1. 実行範囲と結果

job `jobs/nil22_stable_three_lines_t03_actual_curve.json` は実行直前の gate を PASS。
Python 標準ライブラリ Fraction の独立した sparse Laurent 実装により、
全 epsilon と全 m0,…,m5 を同時に保持して ONE batch を実行した。
producer/helper の import、新点、rank/minor/operator、profile、jet は使わない。
制限は 60 秒・1 GiB、実測 0.228746907 秒・最大 RSS 14,592 KiB。
実行 source 全 bytes と SHA は report に埋め込んだ。

確認した全恒等式は次の通り。

- 公開 raw T03 全四 period と H-Casimir gauge の実正逆対応、全 Euler 差。
- source の physical quadric F=k(XY−ZW)、k=−1+epsilon(A−B−C)、t=−k の全式。
- actual GL5 frame とその両 inverse、det=729/epsilon²、tensor scale epsilon。
- 元全 150 係数への forward と inverse、全負 epsilon 冪の消失。
- 端点の四 anchor と full six magnetic bracket の全 150 係数一致。
- 六磁気と period derivatives＋三 translation の全逆、6×6 determinant −3。
- 全 strict 十 Jacobi と全 projective 350 係数の零。
- signed Euler–Pfaffian の全 126 係数、変換則、端点実 D の一致。

全二次括弧は epsilon 次数高々二。D は
`−h²(x+3h)(y−3h)(x−y−3h)` へ厳密に特殊化する。
係数 ring は QQ[epsilon,m0,…,m5,h,Q1,Q2,x,y]、frame 計算時だけ epsilon を反転する。
磁気 parameter は一つも反転せず、全零・residue 零の全境界を含む。

## 2. Whole-source 結論と審査範囲

任意に固定した全六磁気に対し k の定数項は −1 なので、
`epsilon*k != 0` は非空稠密曲線 open。
そこでは公開 true T03 全 source の actual PGL 像と非零 tensor scale である。
多項式性と T03 成分の閉性から、元 full six-magnetic stable_three_lines の全点、
全 PGL 飽和、その閉包が T03 に含まれる。
全源から一点を選んだ推測、境界から逆向きの帰属推測、局所一意性仮定は用いない。

root は実行前に全四 period 正逆、physical F、全六括弧、端点 anchor、
六磁気復帰を独立手査読済み。ここではその手査読と本 producer の全変数算術を区別し、
後続の独立 receipt を未作成のものとして保持する。

## 3. 再生

report 内 `executed_source_bytes` を新しい一時 Python ファイルへそのまま抽出し、
入力 hash を照合した後に同 job の制限下で実行できる。
出力 path は exclusive create なので、再生時は既存 report を上書きせず別作業領域を使う。
依存入力全 bytes の SHA は report の `input_sha256` と同封 manifest に保存した。
report の sparse monomial 順は
`(epsilon,m0,m1,m2,m3,m4,m5,h,Q1,Q2,x,y)`、各係数は QQ の文字列である。
