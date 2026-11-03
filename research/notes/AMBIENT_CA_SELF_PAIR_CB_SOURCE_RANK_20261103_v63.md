# CA self-pair族の実CB lift：全29 source微分はimmersion

2026-11-03。実N/CB incidence packetが与えるCB liftで、全source微分をQQ(a)上で計算した。
結果はrank29であり、明示minorにより全a≠0で成立する。
CBの全fiberやcompleted imageの型はこのrankだけから結論しない。

## 実sourceを先に作ってから微分する

CBの有限sourceは

    f=u0*h²+u1*h*x+u2*x²,
    H=Σ16 cubic_coefficients*cubic_monomial,
    g=mu*f,
    Q=canonical(rawCB(f,mu*f,H))

である。19個のfiber座標とratio muの計20座標を持つ。
muと三つのf係数の積を保持した**実際の非線形source map**を先に作り、20方向全てで微分した。
自由な三つのg係数をratioの代わりに加えることはしない。

今回のliftでは

    f=3a*x², mu=0,
    H=−a*x*w*z+x*w*(h+x)/(3a)+(h+x)²*z/(9a²)+z³/3.

従ってratio微分はcanonical(rawCB(0,3a*x²,0))である。
また全三つのmixed二階微分

    ∂mu∂u_i Q=canonical(rawCB(0,binary_monomial_i,0))

も零に落とさず、全150行で確認した。
全9flag方向は、既知CB parabolic16の補空間から取り、実物理Lie作用で加えた。

## exact rankとminor

mu=0でのfixed-fiber19列は定数QQ行列でrank19。
その単元19×19 blockを消去すると、ratio1+flag9の全150×10行列が残る。
QQ(a) row reductionでrank10、選択minorは

    53747712*a^10/15625

である。従って全sourceの微分rankは29で、全a≠0に対して成立する。
全150×29行列、20座標map、mixed微分、消去した19次inverseと10次minorを保存した。

## 範囲

このpacketは、指定した**一つの実CB lift**がimmersionであることを示す。
全CB proper fiberの支持点数・lengthや、CBの一枝性・smoothness・normalityは示さない。
CAで二つの旗を構成できたことを、CBの二つの旗へ転用しない。
CBの全fiberの独立な判定がない限り、CB completed imageの型は未決のままである。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/ambient_CA_self_pair_CB_source_rank_20261103_v63.py --output <new-output-directory>

code: `code/ambient_CA_self_pair_CB_source_rank_20261103_v63.py`。
全係数とreport: `reports/AMBIENT_CA_SELF_PAIR_CB_SOURCE_RANK_20261103_v63/`。
実行約0.6秒、制限60秒・2.2GB。GBなし。
