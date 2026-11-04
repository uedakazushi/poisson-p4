# 実 T_nil の完全一次 Loewy 層：既存 cache のみ

2026-11-04。**実一次関係のrankは35、HFの初二層は(5,5)、T_nil/m二乗の長さは10**と得た。
一回のcache workerは4.76秒・約50MB、保存25,232 QQ slots・127KBで完了した。

全次数minimal 5×849表示の全一次関係を、既存source cacheと
受理済み定数minorだけから評価する。新しいsource/gauge/IFT、jet、点、GBは行わない。
数値結果は同stemのreportへ記録する。総長、元nil生成元、三乗・指数は別問題である。

## 実source値と正確なfunctional

monicは元11節点のnu=t+s二乗/4である。四つのs/tは保存済み二つの57-source二次解から
読み、CB own-F9を既存full quadratic adapterでCA-F9へ移してell=0を課す。
共有physical9の全405係数、両monicの全一次値とB11像を照合する。
新しい未知数を解く操作は無い。

最終constant組合せcはX0c=0かつH(mXc)0=0を満たす。
残る三wedgeはsector0の(i,7), i=0,4,5のliteral座標であり、canonical m第7行について
m0_7=(5/2)e20、m1_7=r*(保存済みC44の44係数)となる。
従ってwedge一次値の元入力は

  v7,2=(5/2)(Xc)20,2+r*m1_r*(Xc)1

である。ordered-prefix Hadamardを適用し、18pivotの一次消去係数a1に対して
実h_pivot,0*a1を引く。h_pivot,0の対象成分は(5/2)*X_pivot,row20,ziで、未計算商は含まない。
上二targetも同じ実18pivotの一次消去を使う。
最後の81単位Schurでは残る五行のpivot-originが零なので、最終一次値にはこの定数cだけが必要。

## 未知二次値を零に置かない

最初640列のrow20ではResAT0=0、最終全681 cでBold0c=0なので、未知Bold2とadapter2は
寄与しない。既存ResAT1*Bold1だけを保持する。
native140列のrowhの未知proper二次値も、C4偶射影の証明と全681 u5_even2感度零で消える。
従ってraw source二次配列はsurrogateであり、その各列をactual X2と主張しない。
exact constant cで合成した681個の値だけがactualである。

全681本についてm0(Xc)1=0を192係数ずつ検証し、全168 d3一次関係も加える。
得た全849列のrankは固定minorと全列factorizationの両方で証明する。
HF1=40−rank、length(T_nil/m二乗)=5+HF1が正確な結論である。
未知m1を落とした計算でも、d3部分だけの上界でもない。

保存値は全405 shared physical係数、全681×192 source cancellation、固定35minorと
全849列のfactorizationを含む。総長T_nil=L+5の16〜22という別の範囲は維持し、
今回の初二層から正確なLや元nilidealの全生成元・三乗・指数を推論しない。
