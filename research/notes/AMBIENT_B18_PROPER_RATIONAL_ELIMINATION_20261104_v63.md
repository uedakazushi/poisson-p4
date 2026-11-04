# 十四変数 source の五 unit 有理消去と、十二の実 flag 商多項式

2026-11-04。一回の exact builder。新 jet、GB、次の消去は行わない。
元の実 E4 proper-source17式を F0,...,F16 とし、係数順は凍結 packet に固定する。

## 1. 五つの明示解と全 quotient の等号

D=1-d²、B=d(a+d-v2)/2 と書き

    v0=2d-a,
    N=B²+d²B/2-a*d²(d-v2)-d²v4+d*v6,
    e=N/D,
    v1=B/2+2e,
    b=B,
    c=e(a+d+v2)/2-d*v3

とする。新しい変数分母は D だけで D(0)=1。定数除法は2の冪であり、
元の式の既知 QQ 係数（分母3,5等）もそのまま保持する。

単に解を代入して零を確認するだけでなく、順序

    F2 : v0、F0 : b、F4 : v1、F3 : e、F1 : c

で、先行 pivot を消した式が、それぞれ新しい pivot−root の

    1/3, 2, 2, -D, 2

倍であることを厳密に検算する。F4 を使う段階では e をまだ変数として残す。
よって D と非零 QQ 定数だけを反転した環で、選択五式の ideal はこれらの解の graph ideal に等しい。
他の未知 determinant を反転したり、generic な解だけを使ったりしない。

## 2. 十二の実 source 恒等式と flag 因子

残る九変数は (a,d,v2,v3,v4,v5,v6,v7,v8)。
元の残り十二式 F5,...,F16 について、各 term の D 分母次数を数え、
大きな積を展開する前に numerator term 数の保守的上界を保存して150000以下を確認する。
各式の既知 D 冪と QQ 分母を払った numerator を実際に構成する。
その全てが polynomial として正確に d で割れること、および元の有理式へ戻る恒等式を検算する。
primitive な十二整数多項式を H0,...,H11 とすれば、各々

    Fi(root substitution)=constant_unit*d*Hi/D^mi

となる。全 constant と mi を別添に固定する。
従って実 normalized tensor の全完成環は正確に

    QQ[[a,d,v2,...,v8]]/(d*H0,...,d*H11)

である。これは E4 であり actual D0 ではない。
既に証明された proper-source gauge 同値と、上記の明示的 unit quotient を用いる。

## 3. 次の五 pivot は検査のみ

十二 Hi の原点12x9 Jacobian、rank、独立 minor と pivot 分布を報告する。
flag a,d と残る CAfiber のどれを含むかも保存する。
この段階では次の五 pivot を消去しない。
多項式 ideal の GB・local standard basis・全 syzygy・CM・socle・深さは計算しない。
既存四 core への全座標輸送や graph 冗長性の証明も追加しない。

再現:

    .venv/bin/python code/ambient_B18_proper_rational_elimination_20261104_v63.py --output work/B18-proper-rational-replay

全 expanded term150000、出力32MiB、実行60秒・2GBを上限とする。
