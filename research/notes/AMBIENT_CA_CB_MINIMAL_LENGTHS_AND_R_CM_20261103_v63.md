# CA/CB実交叉のminimal成分の一般長とR上のCM・正規な点

2026年11月3日。Dを被約な全射影成分CA、CBのscheme-theoretic交叉とする。
退化五次因子を意味するDではない。既知の被約支持はN27∪R19であり、
今回の全proper模型により(R∩N)red=B18と決定されている。

## 1. 二つのminimal成分の一般長

DのNに沿う一般長は1、Rに沿う一般長は2である。

Nについて、共通一次像Zの一般点の函数体へ平坦に底変換する。
既証の全rank3 first-jet fiber定理は、被約CAとCBの実ファイバーを二つの埋込A10、
その交叉をA9として与える。これは一般Z上の全ファイバーの線形表示であり、
一つの点でのtangentの比較ではない。函数体への平坦な底変換は二つのidealの和と可換し、
実交叉の一般ファイバーも被約A9となる。NはZを優越するので、DのN一般点での長さは1。

Rについては、既証の全R\\Nでの完成環

    Dhat=C[[u1,...,u19,epsilon]]/(epsilon²)

を用いる。その稠密開でnilradicalはepsilon、一般横断長は2である。
従って二つのminimal成分の一般長は正確に(1,2)となる。
埋込素idealやB18でのprimary構造は、この二つの数だけでは決まらない。
これを通常のproper intersectionのChow積として同一視もしない。

## 2. R上のCohen–Macaulayな点は正確にR\\B18

R\\B18=R\\Nの全閉点でDは上のordinary doubleなのでCM、
その被約化Dredは滑らかなRと一致する。

B18の任意の閉点では、DとDredの両方にminimal成分NとRが残り、
閉点でのそれぞれの局所次元は27と19で異なる。
優秀な有限型局所環の完成でもその次元を持つminimal成分が残る。
完全局所CM環はequidimensionalなので、DもDredもそこでCMではない。
従ってR上のCM locusは、DとDredの両方で正確にR\\B18である。

DredのR上のnormal locusも正確にR\\B18である。
開側は滑らかであり、B18側の局所環は相異なる二つのminimal成分を持つので、
normal局所環にはなれない。
一方D自体はR\\B18で非被約、B18では二つのminimal成分を持つため、
R上にnormalな点を一つも持たない。

これはDのR上の全閉点についての結論である。N\\Rの全CM locusや、
R自身のCM locusをこの記述から求めたものではない。
特にDredがB18で非CMであることと、R単独の局所環の性質は区別する。

## 依存証明

- `S61_CA_CB_FULL_RANK3_FIBERS`：全一般Z上の埋込A10、交叉A9。
- `R19_ENTIRE_COMPLEMENT_ORBITS_AND_ORDINARY_DOUBLE`：全R\\Nの実ordinary double。
- `CA_CB_GLOBAL_REDUCED_INTERSECTION_N27_R19`：全minimal支持の尽くし。
- `AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63.md`：全被約交叉R∩N=B18。

新しい大規模計算は不要であり、入力証明と独立監査のhashを
`reports/AMBIENT_CA_CB_MINIMAL_LENGTHS_AND_R_CM_20261103_v63.json`に保存する。
