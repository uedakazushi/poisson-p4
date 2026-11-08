# P*の有限MPS監査：bounded rational quadratic-jet foliationのproper条件

2026-11-08。PRIMARY E_multiple_hyperplane。次の有限literal全source batchに先行する純幾何。
新しいPoisson点／周囲operator／rank／normal deformation jetは用いない。

V4とd≥1に対し、cubic horizontal one-form alpha とGr(2,Sym^d V*)のincidence

```
alpha wedge (G*dF-F*dG)=0
```

はclosed。F,G独立ならGdF-FdGはchar零で非零であり、J=F/Gは非定数rational first integral。
Grassmannianはproperなので、このincidenceのalpha側imageはclosed cone（alpha零も含む）。
共通factor付きpencilも含み、explicit polynomial F/Gがboundaryでdependentになっても
Grass上の二次元spanのproper limitを取るので失われない。

order二以上zero qにおけるquadratic leading bivectorのcubicone-form alpha_Qへ同条件を課し、
properなq∈P4を忘れるとwhole projective tensor space上のclosed条件になる。
従って全sourceの稠密openでdegree固定dの実rational first integralを持つ場合、
sourceの全component closureはこのclosed条件を満たす。

P*は唯一order二zero pのalphaが非零で、非定数rational first integralを持たないことが
`UNIMODULAR_POLE_GLOBAL_MARKED_JET_AND_T36_T38_EXCLUSION.md` §§4–5で証明済み。
よってこの固定d条件をwholeに持つ成分はP*を含めない。

有限batch対象は次の既sourceと幾何から固定したcoordinate zeroだけ：
T24:q=x2、T26/T27/T31/T32/T34/T35:q=x0、T33:q=x3。
全true completed sourceの全parameterを保持し、そのpointの0/1jet零をliteralに確認し、
二次leading Q、cubic alpha、PfQを再構成する。
Pf因子Fiに対しalpha=sum lambda_i*(Pf/Fi)dFiという全coefficient恒等式が得られ、
ratio lambda_i/lambda_jが固定rationalなら、分母を払った整数exponentからF/Gを構成して
上のwedge恒等式と非定数性を直接確認する。
source denominator又はgeneric因子chartだけを反転し、其の全補集合はproper条件により戻す。

全parameter symbolicのこの恒等式は、数字の一点をgenericと呼ぶ方法とは異なる。
Pf零、log表示不成立、residue比非rational又は非一意の場合はそのcaseを未解決と保存し、
追加point／field extension／profile／M1／jetへadaptiveに拡張しない。
本noteの時点では八sourceの陽性又は陰性の算術結論をまだ宣言しない。
