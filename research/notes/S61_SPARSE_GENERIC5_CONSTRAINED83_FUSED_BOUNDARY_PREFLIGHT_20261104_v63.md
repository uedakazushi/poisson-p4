# 29段後の logical final と scratch を分ける読取り専用 preflight

2026-11-04。入力は凍結済み207×54、143215項の残余だけである。
新しい pivot、多項式係数の積・差、gcd、級数計算は行わない。
旧467個の unsafe forecast と29段の算術結果を変更しない。

## 区別する四つの上界

各候補について coefficient を作らず、積 support と和集合を全更新で数える。
一成分の old項数をa、積 support上界をk、結果 support和集合をbとする。
行列そのものの係数上界Lは b−a で更新し、factor列を消去した際にその項数を引く。
最後に退役pivot行R項を引いたものが FINAL logical matrix の上界である。
これは相殺を使わない上界であり、上界超過から実際の結果項数超過を推測しない。

従来の polynomial product/add/back-add コピーは L+2k+3b+64 であった。
これと、matrixだけの途中最大、final matrixだけの上界を別々に保存する。
読取り自体の support-set scratch も、実際に常駐する元143215項と合わせて
150000以内か確認し、時間45秒またはこの条件で止める。

## 本当に fused な exact sparse update の候補

次の方式は数学・予算だけの評価で、未実行である。

1. old成分の係数を一つの monomial→QQ ordered map に入れる。
2. factorと正規化pivot成分の term pairs を直接そのmapへ加減する。
   product多項式、shiftを揃えた全コピー、別の差多項式を同時保持しない。
3. productを一度足し戻してmapをold成分と全係数比較し、再び引いてnewへ戻す。
   これで back-add を丸ごと別多項式として保存せず厳密な恒等式を確認できる。
4. mapから一項ずつ取り出し、FLINT出力へ移してそのmap項を解放する。

mapの項数は常にsupport和集合b以下である。入力old、factor、pivot成分は既存の
行列L内にある。mapと変換中出力は合計b以下に保ち、少量のQQ accumulatorを96項
reserveに含めれば L+b+96 が coefficient/scratch 上界となる。
scalar QQの exact arithmetic と符号付きLaurent exponentを使い、数値代入や
非単項式除算は不要である。これはハッシュを係数恒等式の代用にする案ではない。

## pivot行を外部化する別 representation

必要なら正規化pivot行をliteral receiptへ保存して行列から解放し、一成分だけを
読みながら更新できる。次の二つの走査を、fused版に限って別々に予測する。

- row-major：pivot entryをtarget行ごとに読み、factor列は各行の最後に消す。
- column-major：pivot entryを列ごとに一度読み、factor列は全列の最後に消す。

両者とも各new成分は同じold成分、変更されないfactor列、同じ正規化pivot行から
一度だけ作られるので、走査順を変えても同じ厳密な行操作である。
常駐pivot行R項の代わりに、一成分のnp項だけを数えた全状態上界を保存する。
外部化は未知の値を捨てる操作ではなく、完全なliteral値とsource参照を保つ。
新しいRAW-byte forecast、IOと60秒予算、全係数確認の実装が別に必要で、
このsource-only結果だけでは worker を実行したとはしない。

## common content の構造的 screen

Laurent座標冪を払ったsliceに非零parameter単項式が含まれれば、その全成分の
共通非単項式因子は存在しない。従って各eligible pivotの行・列には、そのpivot
自身があるだけで非単項式contentの可能性が消える。単項式contentを外しても
項数は減らない。

単項式を含まない非零行・列だけを potential slice として列挙する。これは
非自明因子の発見や除算証明ではない。gcdやストリーミング除算の係数増大上界は
未証明であり、新しいgcd jobを起動しない。

形式級数方式も引き続き preflight のみである。actual entry を specialized零で
省略せず、全inverse denominatorと精度伝搬を証明する必要がある。本ノートでは
新しい級数・非単項式pivot・上限変更を実行せず、研究全体の障害とも呼ばない。
