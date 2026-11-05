# 固定 primitive B⁴ 商：四回 lex50ad 長除算の未実行 custody

2026-11-05。ROOT独立五jetで正確なvB=4が受理された後の準備。
canonical4811整数と原D197をheader変更だけで用いる。今回の実行はcompileと`--help`のみで、商値は未計算。

\(Q_2=Q_1/B^4=T_1/(a^3C^4B^4)\) を整数numerator/原D197で得る固定方式。
このQ2は元のS61二次変数q2ではない。Bはprimitive（9−8=1）なのでGaussが整数性を保証する。

## 四回の固定六項除算

lex a>e>dの先頭は50ad。先頭a,d指数を確認し、整数係数を50でexact除去・逆乗算した後、
固定六項の全積をremainderから引く。各回の全remainder零を要求する。
元の4811vectorを保持せず、入力を直接Rへ読み、R/Qの二mapだけで進む。
完成Qをswapして次のRにし、コピーを作らない。

段階支持上限10013→8400→6970→5712→4615、全stageL1境界はSOURCEの有理数式と文字通り同じ。
途中R上限363 bit、最後の全B⁴ multiplyback上限386 bit、hard518。
main52826、最終検証45706 slotsで、全numeric metadata/scratch16000を含む。
既知Bの先頭LC50だけを使い、未知多項式b1を反転しない。

四回の後、固定B⁴をcurrent/next factorだけで作り、old factorを解放してから、
同じimmutable4811入力をRへstream再読する。全係数 \(B^4J_2=I_1\) を検証する。
再読は同じheader/D/order/L1、telemetryで確認し、driverが全SHA256 pinsを実行前後に照合する。
oracle照会や旧572 sourceの再実行ではない。

中間多項式をファイルへ重複保存せず、生成順のFNV64 digestを小さいtelemetryとしてのみ記録する。
これは証明hashではなく、数学証明は各余り零と最終全係数恒等式である。
metadata hash処理は係数product上限1249351と区別し、全metadata領域とCPU制限に含める。

## 保存と停止

`Q1B4QP001`は全current R/Q、fixed/next factor、D/L/Lrとactive quotient、
pass・次のB項・読込位置・全cursor/counterを保持する。
active state1はQ挿入前、2は挿入後のB項更新中なので、半分のremainder状態を曖昧にしない。
full prefix上限5950210 bytes→6 MB、全準備history/入力/code≤1.8 MB、receipt0.2 MBでRAW8 MB。
45秒soft・60秒hard・2 GB。STOPは完成証明ではなく自動retryをしない。
I/O failure時に二つ目の全prefixを作らない。

実際の商計算には新しいROOT source/code/whole auditとONE GOが必要。
他因子・GCD・tail、互いに素、実bad Fitting schemeや補集合消滅、全FULLU57/239は主張しない。
