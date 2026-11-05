# Cprime/ell³ の既知三次 initial closure：独立 SOURCE 監査

2026-11-05。三次の値・rank は計算せず、元 source・monomial quotient・未実行 worker と保存量を検証する。

実際の二次作用素の kernel 列18の像は `(0,275/96,35/48,13/300)` である。その full 二次 initial の pure-plane 部分 P を、保存された全四十四列の九係数と同列の実際の定数 source cofactor から再構成する。P は六項あり、単なる三項の標準 NF representative と混同しない。NF4 は nonplane 成分を読まず、kernel は既知一次式の一次倍の全131次元空間なので、P と full initial の差はその空間に入る。従って P 自身も本物の二次 initial 関係である。

NF4 の像ベクトルを消し、三本の NF3 を作る。標準二次座標は `[44,133,134]`、係数比は `14/55,104/6875` である。全396既知二次関係と P の消滅、単位 embedding を独立に検査する。ann 八十四列は既に旧 NF4 の kernel にあるので、その三次倍756列は既知1980列に含まれる。

495本の三次単項式を、ソートした最初の二変数の二次 divisor と残りの変数に分ける。各 generator について `Id-EQ` が `S1*ker(NF3)` に入る恒等式を普通の単項式積で検査する。embedding E は同じ三次単項式へ重なる場合があるため、全27本の `I-QE` consistency 関係を残す。従って1980本の一次式×S2、P×S1の9本、三 target の ell³、27 consistency の全2019列を27行へ送った cokernel は、Cprime/ell³ の既知三次 initial quotient と正確に一致する。

未実行 worker は二 pass で source を再生成し、選択原列と normalized pivot を高々27本だけ保持する。後者を解放してから固定小 minor/inverse を作り、全列の factorization と source hash を検証する。import は独立に `--help` だけで確認し、main は実行しない。古い custody の重複計数は保存し、新しい disjoint envelope `6000+6000+10000=22000<24000` を使用する。

全投影係数の common denominator は13365000000である。source column の l1 norm と射影係数の最大値から、独立な integer bound184057500000を得る。27次以下の minor は1104bit、保守的な派生 reduced Fraction boundは1222bitで、workerの2048bit上限内に入る。三次2019列の rank を求めずにこの評価を行う。

**正の判定だけ**を許す。rank27なら Nakayama により `m³(Cprime/ell³)=0`、実 HF333から長さは高々9となる。既知DVR rank2・torsion指数 `{1,e>=2}` より長さは `7+min(3,e)>=9` なので、長さ9・e=2が従う。既知三次商が残っても、実際の三次関係を尽くしていないため e>2とは結論しない。root の別 GO 前に三次値 worker を実行しない。
