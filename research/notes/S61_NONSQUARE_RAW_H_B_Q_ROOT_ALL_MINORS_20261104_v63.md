# Q共通根上の元の22minor全部のexact確認

Q成分で得た共通根をd=−N/Dと書く。
N=20(900a+324e²+171e−49)、D=27(3e−2)(6e+1)(e+1)。
この式を保存済みfunction-field共通根の全係数と照合する。

各元minor Piのd次数m≤7についてD^m Pi(a,e,−N/D)を多項式として作り、
モニックa2次Qで割る。全22個のliteral quotient/remainderを保存し、別の
整数分母払い演算で入力評価とQ-membershipを全係数検査する。Hにも同じ検査を
行う。分母を使ったfield cancellationを何度も繰り返す方式ではない。

U=adFBCJも同じ根に代入し、remainderとその2次体ノルムを保存する。
Norm(u+va)=u²−q1uv+q0v²は実際のconjugate積とQ-membershipで検証する。
非零ならgeneric根でU≠0であり、Normの零点は有限例外として保持する。

22個が全て0なら、a/D/NormU非零のQの密開集合のgraphは実際のFULLU内B-rankbad
曲線になる。これはB35の元の最大minorへ既知単元で戻すことで分かる。
full cubic rankやfull quartic rankが落ちるとはそこから推論しない。

どれかが非零なら、その非零ノルムの零点と既存分母例外だけがQ上に残る。
L上generic根がU境界にある証明および前段の有限残余と合わせ、B-bad全体は有限
又は空となる。この場合も空性・全域rank3は別問題である。
