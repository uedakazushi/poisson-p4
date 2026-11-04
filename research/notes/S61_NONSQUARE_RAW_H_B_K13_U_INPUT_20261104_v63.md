# Ka13上のFULLU係数とH4飽和の入力

2026-11-04。source49項のFULLU(a,t,d)をd次数3以下の四つのKa13係数へするだけの入力計算。H4とのgcd、degree4 Euclid、matrix Gaussian、B-minorやNorm4は実行しない。

整数t-powersと共通numeric分母DU*tDen^7を使い、raw a次数89以下の四多項式を作る。g13のpseudo-division商・剰余・最高係数冪を保存して、別Python整数convolutionで全49項と全QQ coefficient normalizationを確認する。52QQ slotのRAW上界719472 bytesは1MiB未満で、source-boundの条件付き承認範囲に入る。

次の候補はH4からU-supported primary factorsを除くことである。Ka13が有限étaleであるgateはg13と導関数のQQ Bézoutまたは受理済みfactor receiptへ明示的に束縛する。field branchごとのH4長は4なので、hU=H4/gcd(H4,U^4)はU非零側の因子を元の重複度のまま残す。Ka13のétale gateを用いなければQQ長52の上界によるU^52が安全な別候補で、今回どちらも計算しない。

future Euclidで必要な係数逆元は各p(a)についてxgcd(p,g13)=1を証明してから使う。非単元pivotなら両proper componentを保持し、非零係数を無条件に反転しない。源のdegree4/12、13成分vector、係数heightを入力として保存し、その後のbudget判断へ渡す。

既存のB21係数から、input supportの小さいrow tripleも数えるが、determinantやNorm4は計算しない。22個の展開minorを再び入力へ戻す方法ではない。

## 確定結果とより小さいnamed-root recipe

U係数のnative計算は0.0345秒で完了した。Uは非零で、39QQ項、numerator/denominator最大1451/1409bit。後段の独立Python mul helperでoperand名をloop変数に再使用した誤りがあり、そのdriver・failure receiptを保存して修正した。完成済みnative証拠のSHAを固定して検証だけを0.1001秒で再開し、全49項・四除算・正規化を受理した。native計算の再実行はしていない。

元の49係数を直接比較し、FULLU=a*d*F*B*C*J（QQ単元は1）を確認した。従ってa,F,C,b1=50a+18t²−9t−2を単元と証明できれば、d側の根は0、rB=−(60t−40)/b1、rJ=(6t+4)/(5a)だけになる。H4をこれらの根で評価し、零ならmonic synthetic divisionを重複度だけ繰り返す方法を次のsource recipeとして保存した。非零評価も単元証明を要求し、proper base componentを落とさない。今回factorや根のKa13評価は実行していない。

代替の固定B-minor候補はrow(2,3,5)、九entry計264QQ項である。実d次数行列からminorはd次数5以下（9ではない）。source integer heightだけによるNorm4の13QQ出力RAW上界は924729 bytesと得た。minor/Norm4自体は未計算で、別予算を要する。
