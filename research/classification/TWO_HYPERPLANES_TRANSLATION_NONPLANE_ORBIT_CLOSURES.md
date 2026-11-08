# 非平面rank二平行移動二全源は二つの22軌道閉包

2026-11-08、PRIMARY `E_reduced_nonnormal`。
全source有限化とcanonicalモデルの全係数算術はproducerで再生済み。
この軌道閉包橋と全ambient算術の独立reviewは未完であり、ambient成分とは数えない。

## 1. 固定cofactorの全source次元は九

`TWO_HYPERPLANES_CUBIC_TRANSLATION_RANK_TWO_NONPLANE_SOURCES.md` の全normalizerと必要十分なJacobi消去を使う。
V=x∂u+h∂vをexactに固定した21係数sourceで、
非平面項は(Ah−Bx)(hu−xv)、(A,B)≠(0,0)である。

A≠0,B=0はAの非零尺度と完全八parameter源なので次元9。
A,B両非零は二つの非零尺度と完全七parameter源なので次元9。
A=0,B≠0はH/x交換で前者へ移り、同じ次元である。
各normalizerの逆は表示された非零尺度だけを使うので、これらは単なるパラメータ上界でなく
全supportを覆う有限chartである。A=B=0の平面sourceはここへ含めない。

## 2. 射影cofactorと二超平面flagの全baseは十三次元

二つの超平面因子を順序付きで選び、
projective bivectorのHamiltonian field `X=P#dlog(H2/H1)` の非零射影classを取る。
rank四のsymplectic openではこの場は非零であり、当該型ではtrace零代表Nが

```
rank N=2, N²=0
```

を満たす。二hyperplane covectorsはP(im N^dual)=P1の相異なる二点である。
順序の交換は有限二択で次元を変えない。

五次元V上のこの[N]は、I=imN（二次元）、K=kerN（三次元）、I⊂Kと
projective isomorphism V/K→Iで完全に決まる。
Iの選択はGr(2,5)の6次元、K/Iの選択はP(V/I)の2次元、
射影isomorphismはP(GL2)の3次元なので全[N]は11次元である。
二hyperplane covectorsの相異なるordered pairがさらに2次元を持ち、全baseは13次元。
このflag空間はGL5が推移的に作用する正常形空間でもある。

その一つのflagをcanonical h,x,[V]へ固定すると、tensorの非零尺度を一意に選んでVをexactに固定できる。
h-Casimir二次代表にEuler ambiguityは無い。
実際E∧Lのh行が零なら `hL_j−x_jL_h=0`、LはlinearなのでL=cEとなりE∧L=0である。
従ってこの正規化したprojective fibreは§1の全21係数源の非平面Jacobi支持そのもので、次元は高々9。

以上から、各非平面全sourceのPGL飽和像Sは

```
dim S <= 13+9=22
```

となる。H/x交換や全parameter境界を落とした固定点の次元計数ではない。

## 3. 実22軌道下界と全sourceへの伝播

`TWO_HYPERPLANES_TRANSLATION_NONPLANE_CANONICAL_MODELS.md` の同じ二点kappa=1で、
全25 Euler列・全25 GL5列・scalar列の行列は両方rank47である。
Eulerの実鉛直rankは24なので、projective orbitは47−24−1=22次元。
全八又は七source微分を追加した行列もrank47で、全350微分との積は零である。
非零minorと全疎行列は
`CLASSIFICATION_TWO_HYPERPLANES_TRANSLATION_NONPLANE_COMPLETION_20261108.json` に保持した。

各normalized sourceは既約で、PGL5も連結・既約だから、その飽和像閉包Sも既約である。
Sはこの22次元orbit Oを含み、§2の上界22を持つ。
従って厳密に

```
S=closure(O).
```

これが全sourceのpositive membershipを単一canonical orbitへ縮約するための必要な橋である。
単なる境界点の所属から未知sourceへ逆向きに推論していない。
任意のPGL安定閉集合Cに対し、canonical点がCに入れば全SもCに入る。
canonical点がCの外なら稠密orbit Oも全て外だが、Sの境界まで外とは主張しない。

## 4. 全ambient判定の正確な停止状態

両canonical点で全Dはrank96、実sourceはrank47で、真法線は七つである。
既存の二源finite contractはこの場合の二次mapを実行しないと定めており、その通り停止した。
射影ambient tangent上界は29で、sourceの22軌道から成分性又は未知親の次元は決まらない。
新法線jet・別kappa・別anchorは追加しない。

次の最弱の全源認識は、同じ二canonical点をfaithful S61等で既知閉成分に入れることである。
陽性なら§3で全sourceが伝播する。
全cokernelや局所一意性を解くことは、そのpositive parent membershipに不要である。
新成分登録は別の全ambient証明と全prior監査が得られた場合だけ行う。
