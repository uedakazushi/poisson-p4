# Ka13上のmonicH4係数だけの準備と次のscalar/B21予算

2026-11-04。承認されたH4係数だけをsource boundに計算する。m2/m3 scalar値、B21 normal form、4×88全matrix、Gaussianは計算しない。

e=tNum/tDenとa^-1=-(g13のa-tail)/g13_constantを使う。Hの74項のnumeric QQ分母をDHで揃え、共通分母DH*g0^5*tDen^8で五つのd係数を整数a多項式にする。各係数をg13でpseudo-divisionし、raw numerator、全商、剰余、最高係数冪を保存する。最後にnumeric contentをcanonical化して五つのKa13係数を得る。別Python整数convolutionで全source substitutionと全除算・正規化恒等式を再検査し、d4係数1を確認する。未知a多項式の逆元は使わない。

次のm2/m3については、既存G商のe次数21/23とa次数、source QQ係数sum、実t/g13 heightから、代入raw次数・pseudo-reduction段数・13QQ出力・全reduction witness保存量の上界だけを導く。実際のscalar値は評価しない。

prepared integer B7×3の21成分は345項、a/e/d次数3/4/3である。d次数3以下なので、そのKa13係数への代入にはH4 reduction自体が不要。既存のknown FULLU row/column-unit scalingを保持したまま、最大1092QQ slotのinput costを別に記録する。22個のmaximal minor2963項を先に展開する方法と混同しない。B21自体のnormal formも今回未計算である。

一回60秒/2GB/32MiB/150k項。全13共役を含む代数上のsource計算であり、新しい数値core点のprobeではない。G0やactual badの実在・空性の結論はここに含めない。

## 確定結果と次の予算候補

承認されたH4-only計算はnative0.1135秒、全体1.8029秒で完了した。四つの非leading係数は計52個のQQ係数で、monic leading1を含めると53項。実numerator最大1986bit、denominator最大1922bitである。74 source全項の代入、五つのpseudo-division、canonical normalizationを別の整数演算で全係数照合した。

m2/m3については評価せず、全source/reduction certificateを含む保守的RAW上界をそれぞれ8836406/10542193 bytes、合計19378599 bytesと得た。13QQ出力そのものは各561886/616603 bytes以下である。従って次に一回の二scalar評価を提案できるが、このpacketはその実行許可ではない。

prepared B21のKa13係数は最大1092QQ slot、numerator/denominator上界12334/12008bit、RAW8165976 bytes以下。d次数3のためここにはH4還元は不要である。この別入力を先に作る選択肢も保存し、全4×88 matrix/Gaussianへ自動で進まない。
