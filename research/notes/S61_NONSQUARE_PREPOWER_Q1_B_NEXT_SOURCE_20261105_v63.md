# B付値4の一次結果、独立五jetと固定B⁴商のSOURCE

2026-11-05。一次report `S61_NONSQUARE_PREPOWER_Q1_B_HASSE_ENTRY_20261105_v63`
SHA `b573db2d9d165cd311786faadf3f2dd1ca3488276604086bdd15e396e7d509a7`。
k0..3全零、k4は1107項非零で、一次結果は正確な \(v_B(Q_1)=4\)。
全三pass/readback、cache3604、peak35044/37436、全約661万products、RAW1123738 bytes、約1.02秒。
まだ独立五jetの値検証とB商は実行していない。

## 別算術：五jet下降d-Horner

元4811をd次でI_lに分け、H0=I18,H1..4=0から、l17..0、r4..1降順で
\[
H_r\leftarrow R H_r+b_1H_{r-1},\qquad
H_0\leftarrow R H_0+b_1^{18-l}I_l
\]
を行う。全五層の分母はD b1¹⁸。binomial/S-cache/一次nativeを使わず、
四層全零と最初層1107全係数を照合するSOURCEである。
五jet合計11295、一temp2297、current/next累乗685、input4811とmetadata16000込み54176 slots。
L1上限436、preop437、hard518、上限300万products、RAW8 MB。
b1を反転せず、b1=0のsource成分を捨てない。

## 固定六項Bによる商：値未計算

Bは整数primitiveなので、独立受理後には \(J_2=I_1/B^4\) は整数多項式で原D197を保てる。
lex a>e>dの先頭単項式は50ad。固定LC50の整数長除算を四回行い、
先頭a,d指数と係数50可除性、各全余り零を要求できる。未知多項式b1の反転は不要。

真の段階支持上限は10013→8400→6970→5712→4615、次数boxも各段階に保存する。
途中Rは旧input−B×真の商の部分和なので旧boxから出ない。
\(B=50ad+U\)、L1(U)=129を用い、m=旧a次数に対する商L1上限は
\(L/50\sum_{h=0}^{m-1}(129/50)^h\)。全四段階の整数L1境界を有理数でSOURCEに固定する。
段階remainder363 bit、最終の全係数B⁴ multiplyback386 bit、hard518。

全所有量52826 slots、最終検証45706、全product上限1249351。
中間多項式は同時保持せずstream hashのみ。最終商だけを保存し、原inputを読み直して
全係数 \(B^4J_2=I_1\) を検証する方式を提案する。
全prefix上限5950210 bytes→6 MB予約、準備入力/コード/全history≤1.8 MB、receipt0.2 MBでRAW8 MB。

これはSOURCEのみで、B⁴商には独立付値peerと新しいworker/code/whole監査・ONE GOが必要。
既知unitによる入力縮小以外の、GCD・互いに素・bad Fitting等号・補集合消滅・全FULLU57/239は主張しない。
